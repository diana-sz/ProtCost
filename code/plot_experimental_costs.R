library(here)
library(RColorBrewer)

setwd(here())

source("code/uni_colors.R")

colors <- c(uni_blue)
colors <- brewer.pal(6, "Paired")[c(2,3,4,5)]
cex_all <- 1.05

# ATP synthase
exp1 <- data.frame(
  titrated_level = c(
    # Glucose
    0, 0.2, 0.4, 0.6, 0.8, 1, 1.2, 1.5, 2, 3, 4,
    # Succinate
    0.1, 0.5, 0.75, 1, 1.5, 2.5, 3, 3.5, 4
  ),
  mu_rel = c(
    # Glucose
    0.75, 0.83, 0.85, 0.9, 0.925, 0.95, 0.99, 1, 0.95, 0.91, 0.89,
    # Succinate
    0.48, 0.87, 0.95, 1, 0.95, 0.8, 0.7, 0.625, 0.52
  ),
  condition = c(
    rep("Glucose", 11),
    rep("Succinate", 9)
  ),
  experiment = rep(1, 20)
)


# citrate synthase
exp2 <- data.frame(
  titrated_level = c(
    0.0, 0.2, 0.5, 1.0, 2.0,       # Glucose/Acetate
    0.0, 0.1, 0.2, 0.4, 0.8, 1.6, 2.4  # Acetate
  ),
  mu_rel = c(
    0.615384615, 0.953846154, 0.969230769, 1.0, 0.984615385,
    0.428571429, 0.514285714, 0.571428571, 0.771428571, 1.0, 0.971428571, 0.942857143
  ),
  condition = c(rep("Glucose/Acetate", 5), rep("Acetate", 7)),
  experiment =  rep(2, 12)
)


# L. lactis
exp3 <- data.frame(
  titrated_level = c(
    # LDH
    0.0, 0.1, 0.6, 0.9, 1.0, 1.1, 1.2, 1.3,
    # PFK
    1.0, 1.5, 3.6, 4.0, 6.0,
    # LAS
    0.75, 1, 1.5, 2, 4, 6,
    # GAPDH
    0.125, 0.25, 0.75, 1, 1.5, 2, 2.5
  ),
  mu_rel = c(
    # LDH
    0.72, 0.75, 0.94, 0.93, 0.98, 0.95, 0.95, 0.96,
    # PFK
    1.0, 0.986666667, 0.973333333, 0.933333333, 0.933333333,
    # LAS
    0.875, 1, 1.016666667, 0.966666667, 0.916666667, 0.875,
    # GAPDH
    0.533333333, 0.906666667, 0.973333333, 1.013333333, 1.0, 1.026666667, 1.013333333
  ),
  condition = c(
    rep("LDH", 8),
    rep("PFK", 5),
    rep("LAS", 6),
    rep("GAPDH", 7)
  ),
  experiment = rep(3, 26)
  
)


# transporter
exp4 <- data.frame(
  titrated_level = c(
    0,
    0.111111111,
    0.222222222,
    0.444444444,
    0.666666667,
    0.888888889,
    1.333333333,
    2.666666667
  ),
  mu_rel = c(
    0.46875,
    0.625,
    0.8125,
    0.9375,
    0.9375,
    1.0,
    1.03125,
    0.90625
  ),
  condition = rep("Maltose", 8),
  experiment =  rep(4, 8)
  
)

all_data <- rbind(exp1, exp2, exp3, exp4)


glucose <- function(x) 0.75 + 0.7 * (x / (0.85 + x))^(1.2) - 0.11 * x
succinate <- function(x) 0.35 + 1.6 * x / (0.7 + x) - 0.3 * x
glucose_acetate <- function(x) 0.6 + 6.25 * ( (x/0.15)^3 / (1 + (x/0.15)^3) ) * 1 / (15 + x)
acetate <- function(x) 0.43 + 2.7 * ( (x/0.4)^2 / (1 + (x/0.4)^2) ) * 1 / (3 + x)
LDH <- function(x) 0.7 + 0.5 * x / (0.7 + x) * 1 / (1 + (x/1.5)^7)
PFK <- function(x) 1.01 - 0.015 * x
LAS <- function(x) 0.247589 + 11.5045 * x^2 / (0.419085^2 + x^2) / (12.6873 + x)
GAPDH <- function(x) -1.4333 + 602.464 * x^2 / (0.0622921^2 + x^2) / (244.484 + x)
maltose <- function(x) 0.15/0.32 + 1.25 * (x / (0.5 + x))^1.2 - 0.22 * x

x <- seq(0, 4, length.out = 50)
exp_curves <- list(
  "Glucose"         = list(x = seq(0, 4.5, length.out = 50),  y = glucose(seq(0, 4.5, length.out = 50))),
  "Succinate"       = list(x = seq(0, 4.5, length.out = 50),  y = succinate(seq(0, 4.5, length.out = 50))),
  "Glucose/Acetate" = list(x = seq(0, 2.5, length.out = 50),  y = glucose_acetate(seq(0, 2.5, length.out = 50))),
  "Acetate"         = list(x = seq(0, 3,   length.out = 50),  y = acetate(seq(0, 3,   length.out = 50))),
  "LDH"             = list(x = seq(0, 1.4, length.out = 50),  y = LDH(seq(0, 1.4, length.out = 50))),
  "PFK"             = list(x = seq(0.8, 4, length.out = 50),  y = PFK(seq(0.8, 4, length.out = 50))),
  "LAS"             = list(x = seq(0.4, 4, length.out = 50),  y = LAS(seq(0.4, 4, length.out = 50))),
  "GAPDH"           = list(x = seq(0, 3.5, length.out = 50),  y = GAPDH(seq(0, 3.5, length.out = 50))),
  "Maltose"         = list(x = seq(0, 3,   length.out = 50),  y = maltose(seq(0, 3,   length.out = 50)))
)

                   

png(paste0("figures/experimental_costs.png"),
    type="cairo", units="cm", res = 300,
    width=17, height=13)


par(mfrow=c(2,2), oma = c(1,1,0.5,0))

for(exp in unique(all_data$experiment)){
  one_exp <- all_data[all_data$experiment == exp,]
  left <- ifelse (exp %in% c(1,3), 2.5, 1)
  par(mar = c(2.5, left, 1, 3.5-left))
  
  conditions <- unique(one_exp$condition)
  for(cond in seq_along(conditions)){
    one_cond <- one_exp[one_exp$condition == conditions[cond],]
    
    if(cond == 1){
      plot(
        mu_rel ~ titrated_level, data = one_cond, 
        axes = FALSE,
        xlim = c(0, 4), ylim = c(0, cex_all),
        pch = 19,
        col = colors[cond],
        cex = cex_all, cex.lab = cex_all, cex.axis = cex_all
      )
    }else{
      points(one_cond$titrated_level, one_cond$mu_rel, 
            pch = 19,
            col = colors[cond],
            cex = cex_all, )
    }
    
    curve_info <- exp_curves[[conditions[cond]]]
    lines(curve_info$x, curve_info$y, col = colors[cond], lwd = 1.5)
    
    mtext(
      paste0("(", letters[exp], ")"),
      side = 3,
      adj = 0.5,
      line = 0.35,
      cex = cex_all,
      font = 2
    )

    labels <- if (exp %in% c(1,3)) c("0", "", "0.5", "", "1") else FALSE
    axis(1, at = 0:4)
    axis(2, at = c(0, 0.25, 0.5, 0.75, 1), labels = labels)
    box()
  }

  legend("bottomright", legend = conditions, col = colors[1:length(conditions)], pch = 19, bty = "n", cex=cex_all)

}

mtext("Titrated level/wild type level", side = 1, outer = TRUE, line = -0.3, cex = cex_all)
mtext("Relative growth rate", side = 2, outer = TRUE, line = -0.3, cex = cex_all)

dev.off()

