#!/usr/bin/env Rscript
# Real survival analysis via the survival package, called by
# ultimate.survival_backend when lifelines is unavailable. This script writes a
# genuine KM table and figure (via ggsave/base plotting), never a fixed curve.

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default = NA_character_) {
  idx <- which(args == flag)
  if (length(idx) == 0 || idx[1] == length(args)) return(default)
  args[idx[1] + 1]
}

input_path <- get_arg("--input")
output_path <- get_arg("--output")
figure_path <- get_arg("--figure")
label <- get_arg("--label", "feature")

if (is.na(input_path) || !file.exists(input_path)) {
  cat(sprintf('{"status":"skipped","reason":"missing_input"}\n'))
  quit(status = 0)
}
if (!requireNamespace("survival", quietly = TRUE)) {
  cat(sprintf('{"status":"skipped","reason":"dependency_missing:survival"}\n'))
  quit(status = 0)
}

data <- read.delim(input_path, stringsAsFactors = FALSE)
data$time <- suppressWarnings(as.numeric(data$time))
data$event <- suppressWarnings(as.numeric(data$event))
data <- data[!is.na(data$time) & data$time >= 0 & !is.na(data$event), ]
if (nrow(data) < 2 || length(unique(data$group)) < 2) {
  cat(sprintf('{"status":"skipped","reason":"insufficient_samples_for_two_groups"}\n'))
  quit(status = 0)
}

library(survival)
fit <- survfit(Surv(time, event) ~ group, data = data)
km <- summary(fit)
km_table <- data.frame(
  time = km$time,
  group = rep(names(fit$strata), km$strata),
  surv = km$surv,
  lower = km$lower,
  upper = km$upper
)
write.table(km_table, output_path, sep = "\t", row.names = FALSE, quote = FALSE)

lr <- tryCatch(survdiff(Surv(time, event) ~ group, data = data), error = function(e) NULL)
logrank_p <- if (!is.null(lr)) 1 - pchisq(lr$chisq, df = length(lr$n) - 1) else NA_real_

hr <- hr_low <- hr_high <- NA_real_
tryCatch({
  cx <- coxph(Surv(time, event) ~ group, data = data)
  hr <- exp(unname(coef(cx))[1])
  ci <- exp(unname(confint(cx))[1, ])
  hr_low <- ci[1]
  hr_high <- ci[2]
}, error = function(e) NULL)

png(figure_path, width = 1400, height = 1000, res = 150)
plot(fit, xlab = "Time", ylab = "Survival probability",
     main = sprintf("Kaplan-Meier: %s", label), col = c("#2E6F9E", "#B5651D"))
legend("topright", legend = names(fit$strata), col = c("#2E6F9E", "#B5651D"), lty = 1)
if (!is.na(logrank_p)) {
  mtext(sprintf("Log-rank p = %.4g", logrank_p), side = 3, line = 0.4)
}
invisible(dev.off())

cat(sprintf('{"status":"ready","engine":"r_survival","logrank_p":%s,"cox_hr":%s,"cox_hr_low":%s,"cox_hr_high":%s}\n',
            ifelse(is.na(logrank_p), "null", format(logrank_p, digits = 10)),
            ifelse(is.na(hr), "null", format(hr, digits = 10)),
            ifelse(is.na(hr_low), "null", format(hr_low, digits = 10)),
            ifelse(is.na(hr_high), "null", format(hr_high, digits = 10))))
