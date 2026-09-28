# Simulation-based minimum detectable moderation effects
# Each model uses 200 pilot and 2,000 main replicates; alpha = .05 and power = .80.
# Failed fits remain in the power denominator.

simulate_moderation_power <- function (k, obj, od) 
{
    id <- names(obj$models)[k]
    original <- obj$models[[id]]
    target <- obj$targets[[id]]
    d <- model.frame(original)
    ff <- formula(original)
    environment(ff) <- environment()
    m <- lmer(ff, data = d, REML = isREML(original), control = lmerControl(optimizer = "bobyqa", 
        optCtrl = list(maxfun = 2e+05)))
    j <- match(target, names(fixef(m)))
    set.seed(202609100 + k)
    b <- fixef(m)
    b[j] <- 0
    nsim <- 2200L
    sim <- simulate(m, nsim = nsim, newparams = list(beta = b, theta = getME(m, "theta"), sigma = sigma(m)), 
        re.form = NA)
    one <- function(i) {
        warn <- character()
        err <- NA_character_
        mm <- tryCatch(withCallingHandlers(as_lmerModLmerTest(refit(m, sim[[i]])), warning = function(w) {
            warn <<- c(warn, conditionMessage(w))
            invokeRestart("muffleWarning")
        }), error = function(e) {
            err <<- conditionMessage(e)
            NULL
        })
        if (is.null(mm)) 
            return(data.frame(replicate = i, estimate = NA, SE = NA, df = NA, singular = NA, 
                warning = paste(warn, collapse = " | "), error = err))
        co <- coef(summary(mm))[target, ]
        conv <- mm@optinfo$conv$lme4$messages
        conv <- conv[!grepl("boundary.*singular", conv)]
        bad <- length(conv) > 0 || any(grepl("failed to converge|negative eigenvalue|unable to evaluate|degenerate|Hessian", 
            warn, ignore.case = TRUE))
        data.frame(replicate = i, estimate = unname(co["Estimate"]), SE = unname(co["Std. Error"]), 
            df = unname(co["df"]), singular = isSingular(mm), convergence_failed = bad, warning = paste(unique(c(warn, 
                conv)), collapse = " | "), error = NA_character_)
    }
    rr <- vector("list", nsim)
    for (i in seq_len(nsim)) {
        rr[[i]] <- one(i)
        if (i%%200 == 0) {
            saveRDS(rr[seq_len(i)], file.path(od, paste0(id, "_checkpoint.rds")))
            cat(id, i, "/", nsim, "\n")
        }
    }
    rr <- bind_rows(rr)
    write_csv(rr, file.path(od, paste0(id, "_replicates.csv")))
    finite <- with(rr, is.finite(estimate) & is.finite(SE) & is.finite(df) & SE > 0 & df > 
        0 & !is.na(convergence_failed) & !convergence_failed)
    evaluate <- function(effect, index) {
        r <- rr[index, ]
        ok <- finite[index]
        p <- rep(NA_real_, nrow(r))
        p[ok] <- 2 * pt(-abs((r$estimate[ok] + effect)/r$SE[ok]), df = r$df[ok])
        reject <- sum(p < 0.05, na.rm = TRUE)
        ci <- binom.test(reject, nrow(r))$conf.int
        data.frame(effect = effect, n = nrow(r), usable = sum(ok), failed = sum(!ok), singular = sum(r$singular, 
            na.rm = TRUE), rejections = reject, power = reject/nrow(r), lower = ci[1], upper = ci[2])
    }
    pilot <- seq_len(200)
    main <- 201:2200
    se <- sqrt(vcov(m)[j, j])
    grid <- se * c(0, 1, 2, 2.5, 3, 3.5, 4, 5, 6)
    pg <- bind_rows(lapply(grid, evaluate, index = pilot))
    write_csv(pg, file.path(od, paste0(id, "_pilot_grid.csv")))
    high <- max(grid)
    while (evaluate(high, main)$power < 0.8 && high < 100 * se) high <- 2 * high
    threshold <- if (evaluate(high, main)$power < 0.8) 
        NA_real_
    else uniroot(function(e) evaluate(e, main)$power - 0.8, c(0, high), tol = se/10000)$root
    effects <- unique(c(0, if (is.finite(threshold)) threshold * c(0.9, 1, 1.1) else high))
    ans <- bind_rows(lapply(effects, evaluate, index = main)) %>% mutate(model = id, seed = 202609100 + 
        k, target = target, MDES80 = threshold, reported_MDES80 = abs(obj$unit_factors[[id]]) * 
        threshold, reported_effect = abs(obj$unit_factors[[id]]) * effect, n_subjects = nlevels(getME(m, 
        "flist")[[1]]), n_rows = nobs(m), sigma = sigma(m), random_intercept_sd = as.data.frame(VarCorr(m))$sdcor[1], 
        REML = isREML(m))
    write_csv(ans, file.path(od, paste0(id, "_summary.csv")))
    ans
}

run_moderation_sensitivity <- function (fitted_models, output_dir) 
{
    ids <- c(paste0("EP_", c("hrv_baseline_z", "hrv_mvc_z", paste0(EMA_VARIABLES, "_z"))), 
        unlist(lapply(c("effSens", "rewSens"), function(y) paste("DM", y, c("hrv_baseline", 
            "hrv_mvc", paste0(EMA_VARIABLES, "_z")), sep = "_"))))
    obj <- list(models = lapply(fitted_models[ids], `[[`, "fit"), targets = list(), unit_factors = list())
    for (id in ids) {
        terms <- names(lme4::fixef(obj$models[[id]]))
        target <- if (startsWith(id, "EP_")) 
            grep("^force_z:phaseLuteal:", terms, value = TRUE)
        else grep("^phase(Luteal|1):", terms, value = TRUE)
        require_condition(length(target) == 1L, paste("Ambiguous interaction:", id))
        obj$targets[[id]] <- target
        obj$unit_factors[[id]] <- if (startsWith(target, "phase1:")) 
            -2/stats::sd(stats::model.response(stats::model.frame(obj$models[[id]])))
        else 1
    }
    od <- file.path(output_dir, "moderation_sensitivity")
    dir.create(od, recursive = TRUE, showWarnings = FALSE)
    cl <- parallel::makePSOCKcluster(4, outfile = file.path(od, "simulation.log"))
    on.exit(parallel::stopCluster(cl), add = TRUE)
    parallel::clusterEvalQ(cl, {
        suppressPackageStartupMessages({
            library(lmerTest)
            library(dplyr)
            library(readr)
        })
    })
    parallel::clusterExport(cl, "simulate_moderation_power", envir = environment())
    result <- parallel::parLapplyLB(cl, seq_along(obj$models), simulate_moderation_power, obj = obj, 
        od = od)
    write_table(dplyr::bind_rows(result), od, "minimum_detectable_effects")
    invisible(result)
}

