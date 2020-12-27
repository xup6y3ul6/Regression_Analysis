library(tidyverse)
library(GGally)

# plot ====
## for variable
get_histogram <- function(Data, variable){
  ggplot(Data, aes_string(x = variable, fill = "location")) +
    geom_histogram(color = "white") +
    theme_bw() +
    theme(legend.position = "top")
}

get_scatter_plot <- function(Data, variable){
  ggplot(Data, aes_string(x = variable, y = "agree_rate", color = "location")) +
    geom_point() +
    theme_bw() +
    theme(legend.position = "top")
}

custom_cor_color <- function(data, mapping, color = I("black"), sizeRange = c(1, 5), ...) {
  
  # get the x and y data to use the other code
  x <- GGally::eval_data_col(data, mapping$x)
  y <- GGally::eval_data_col(data, mapping$y)
  
  ct <- cor.test(x,y)
  sig <- symnum(
    ct$p.value, corr = FALSE, na = FALSE,
    cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1),
    symbols = c("***", "**", "*", ".", " ")
  )
  
  r <- unname(ct$estimate)
  rt <- format(r, digits=2)[1]
  tt <- as.character(rt)
  
  
  # plot the cor value
  p <- ggally_text(
    label = tt, 
    mapping = aes(),
    xP = 0.5, yP = 0.5, 
    size = 6,
    color=color,
    ...
  ) +
    # add the sig stars
    geom_text(
      aes_string(
        x = 0.8,
        y = 0.8
      ),
      label = sig, 
      color = color,
      ...
    ) + 
    theme(
      panel.background=element_rect(fill="white", color = "black", linetype = "dashed"),
      #panel.background = element_rect(color = "black", linetype = "dashed"),
      panel.grid.minor=element_blank(),
      panel.grid.major=element_blank()
    ) 
  
  corColors <- RColorBrewer::brewer.pal(n = 9, name = "RdBu")
  # RcolorBrewer::display.brewer.pal(n=9, name = "ReBu)
  if (r <= -0.6) {
    corCol <- corColors[2]
  } else if (r <= -0.3) {
    corCol <- corColors[3]
  } else if (r <= -0.145) {
    corCol <- corColors[4]
  } else if (r < 0.145) {
    corCol <- corColors[5]
  } else if (r < 0.3) {
    corCol <- corColors[6]
  } else if (r < 0.6) {
    corCol <- corColors[7]
  } else {
    corCol <- corColors[8]
  } 
  
  p <- p + theme(
    panel.background = element_rect(fill= corCol)
  )
  
  p
}

custom_smooth <- function(data, mapping, ...) {
  ggplot(data = data, mapping = mapping) +
    geom_point(size = 0.5, color = "dimgray") + 
    geom_smooth(method = "loess", size = 0.8, 
                color = "forestgreen", fill = "lightgreen", ...) +
    theme_bw()
}

## for GLM
get_fitted_plot <- function(model, .title, .color){
  df <- data.frame(true_value = model$y,
                   fitted_value = model$fitted.values)
  
  g <- ggplot(df, aes(x = fitted_value, y = true_value)) +
    geom_point(color = .color, shape = 1) +
    geom_abline(slope = 1, intercept = 0, color = "tomato", linetype = "dashed") + 
    coord_cartesian(xlim = c(0.15, 0.5), ylim = c(0.15, 0.5)) +
    labs(title = .title) +
    theme_bw()
  
  return(g)
}

get_eta_plot <- function(model, .title, .color, .xlim){
  df <- data.frame(true_value = model$y,
                   eta_value = model$linear.predictors)
  
  g <- ggplot(df, aes(x = eta_value, y = true_value)) +
    geom_point(color = .color, shape = 1) +
    stat_function(fun = binomial_link, args = list(model = .title),
                  xlim = .xlim, color = "tomato") + 
    labs(title = .title) +
    theme_bw()
  
  return(g)
}

get_rstandard_plot <- function(model, .title, .color, .type){
  df <- data.frame(fitted_value = model$fitted.values,
                   standardized_residual = rstandard(model, type = .type))
  
  g <- ggplot(df, aes(x = fitted_value)) +
    geom_point(aes(y = standardized_residual), color = .color, shape = 1) +
    geom_hline(yintercept = qnorm(c(0.025, 0.5, 0.975), mean = 0, sd = 1), 
               color = c("tomato", "black", "tomato"), 
               linetype = c("dashed", "solid", "dashed")) +
    # scale_shape_identity("Residual", guide = "legend",
    #                      breaks = c(19, 1), labels = c("Pearson", "deviance")) +
    labs(title = .title) +
    theme_bw()
  
  return(g)
}

## for LM
get_residual_plot <- function(model){
  data.frame(fitted_value = model$fitted.values,
             residual =  model$y - model$fitted.values, 
             name = rownames(model$data)) %>% 
    ggplot(aes(x = fitted_value, y = residual)) +
    geom_point() +
    geom_hline(yintercept = 0, color = "tomato", linetype = "dashed") +
    theme_bw()
}

get_boxcox_plot <- function(model, lambda){
  boxcox.sse <- function(lambda, model) {
    X <- model.frame(model)[,-1]
    Y <- model.frame(model)[,1]
    n <- nrow(model.frame(model))
    K2 <- prod(Y)^(1/n)            # (3.36a)
    K1 <- 1/(lambda*K2^(lambda-1)) # (3.36b)
    
    if (lambda == 0) {
      W <- K2*log(Y)        # (3.36)
    } else{
      W <- K1*(Y^lambda-1)  # (3.36)
    }
    # Deviance = Residual Sum of Squares, SSE
    .df <- data.frame(W, X)
    return(deviance(lm(W ~ ., .df))) #有問題！！！！！！
  }
  
  SSE <- sapply(lambda, boxcox.sse, lm_w2)
  .data <- data.frame(lambda, SSE)
  
  ggplot(.data, aes(lambda, SSE)) +
    geom_line() + 
    scale_x_continuous(breaks = seq(min(lambda), max(lambda), 0.5)) +
    theme_bw()
}

# link ====

binomial_link <- function(x, model){
  if (model == "binom_identity") {
    return(x)
  } else if (model == "binom_logit") {
    return(1/(1+exp(-x)))
  } else if (model == "binom_probit") {
    return(pnorm(x))
  } else if (model == "binom_cloglog") {
    return(1-exp(-exp(x)))
  } else {
    print("something wrong ><")
  }
}


# likelihood and AIC ====

binom_loglikelihood <- function(model){
  y= model$data$agree_vote
  N = model$data$valid_vote
  fit = model$fitted.value
  sum(log(dbinom(y, size = N, prob = fit)))
}

get_binom_XIC <- function(model, k = 2, Npar){
  ll <- binom_loglikelihood(model)
  xic <- -2*ll + k*Npar
  return(xic)
}
