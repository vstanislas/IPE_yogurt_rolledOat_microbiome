
'%!in%' <- function(x,y)!('%in%'(x,y))



colorize <- function(x, color) {
  if (knitr::is_latex_output()) {
    sprintf("\\textcolor{%s}{%s}", color, x)
  } else if (knitr::is_html_output()) {
    sprintf("<span style='color: %s;'>%s</span>", color,
            x)
  } else x
}


library(MASS)
library(crayon)

#' Transform variables in a dataframe
#'
#' @param DF a data frame.
#' @param trans a named vector specifying the transformation to apply to each variable, where each name corresponds to each variable name. Transformations available are either: log2 (\item{"log"}, \item{"ln"}), pseudo log2 (\item{"psd_log"}), square root (\item{"sqrt"}), Box-Cox (\item{"bc"}), Box-cox with negatives values (\item{"bcn"}), no transformation (\item{""}, \item{" "}) 
#'
#' @return Returns DF with the transformed variables added as new columns where each new column name is the original variable name appended with the transformation name.
transform_data <- function(DF, trans){
  trans_name <- unique(trans)
  
  for(i in seq(length(trans_name))){
    var_to_trans <- names(trans[trans==trans_name[i]])
    
    if(trans_name[i] %in% c("log", "ln")){
      cat(blue(paste("Use log2 transformation on the",  length(var_to_trans), "following variables:\n")))
      cat(var_to_trans, sep="\n")
      new_var <- log_transform(DF[, var_to_trans])
    }
    
    if(trans_name[i] %in% c("psd_log")){
      cat(blue(paste("Use pseudo log2 transformation on the",  length(var_to_trans), "following variables:\n")))
      cat(var_to_trans, sep="\n")
      new_var <- pseudo_log(DF[, var_to_trans], base=2)
    }
    
    if(trans_name[i] == "sqrt"){
      cat(blue(paste("Use sqrt transformation on the",  length(var_to_trans), "following variables:\n")))
      cat(var_to_trans, sep="\n")
      new_var <- sqrt(DF[, var_to_trans])
    }
    
    if(trans_name[i] == "bc"){
      cat(blue(paste("Use Box-cox transformation on the",  length(var_to_trans), "following variables:\n")))
      cat(var_to_trans, sep="\n")
      new_var <- c()
      for(j in seq(length(var_to_trans))){
        x <<- DF[, var_to_trans[j]]
        b <- boxcox(lm(x ~ 1), plotit=F)
        lambda <- b$x[which.max(b$y)]
        ci_lim <- max(b$y) - 1/2 * qchisq(.95,1)
        ci <- c(head(b$x[b$y > ci_lim], 1), tail(b$x[b$y > ci_lim], 1))
        cat(blue(paste(" lambda ",  var_to_trans[j], ": ", round(lambda, 2), " [", round(ci[1], 2), ";", round(ci[2], 2), "]", "\n", sep="")))
        new_var_j <- bc_tr(x, lambda)
        new_var <- cbind(new_var, new_var_j)
      }
    }
    
    if(trans_name[i] == "bcn"){
      cat(blue(paste("Use Box-cox with negatives allowed transformation on the",  length(var_to_trans), "following variables:\n")))
      cat(var_to_trans, sep="\n")
      new_var <- c()
      for(j in seq(length(var_to_trans))){
        pt <- powerTransform(DF[, var_to_trans[j]], family="bcnPower")
        new_var_j <-  bcnPower(DF[, var_to_trans[j]], lambda=pt$lambda, gamma=pt$gamma)
        new_var <- cbind(new_var, new_var_j)
      }
    }
    
    if(trans_name[i] %in% c("", " ")){
      cat(blue(paste("No transformation on the",  length(var_to_trans), "following variables:\n")))
      cat(var_to_trans, sep="\n")
      new_var <- NULL
    }
    
    if(!(trans_name[i] %in% c("log", "psd_log", "sqrt", "bc", "bcn", "", " "))){
      cat(blue(paste("No transformation,", trans_name[i], "not implemented\n")))
      new_var <- NULL
    }
    
    if(!is.null(new_var)){
      new_var <- data.frame(new_var)
      colnames(new_var) <- paste(var_to_trans, trans_name[i], sep="_")
      DF <- cbind(DF, new_var)
    }
  }
  return(DF=DF)
}


bc_tr <- function(x, lambda){
  if(lambda == 0) y <- log(x)
  else  y <- (x ^ lambda - 1) / lambda
  return(y)
}

LOG <- function(x){ ## from Maaslin2
  y <- replace(x, x == 0, min(x[x > 0])/2)
  return(log2(y))
}

log_transform <- function(DF){
  
  if(is.vector(DF)) DF <- LOG(DF)
  else DF <- apply(DF, 2, LOG)
  return(DF)
}

# from scales::pseudo_log_trans and 
# https://win-vector.com/2012/03/01/modeling-trick-the-signed-pseudo-logarithm/
pseudo_log <- function(x, sigma = 1, base = exp(1)){asinh(x/(2 * sigma))/log(base) }





plot1Var <- function(var, title="", subtitle="", cols=c("#0072B2A0", "#00B259A0"), 
                    axis_angle=F, DATA.=NULL, size_counts=3, printPercent=NULL, print="number_pct"){
  
  pp <- ggplot(DATA., aes(x=.data[[var]]))
  
  class_var <- class(DATA.[, var])[1]
  if(class_var %in% c("integer", "numeric")){
    pp <- pp + geom_histogram(bins = 12, color= cols[1], fill= cols[1], alpha=0.4)
  }
  else if(class_var %in% c("factor", "ordered")){
    pp <- pp + geom_bar(color= cols[2], fill= cols[2])
    
    if(!is.null(printPercent)) if(printPercent) print="number_pct" else print="number" # to keep the old code with printPercent working
    
    if(print=="number_pct"){
      pct_format = scales::percent_format(accuracy = .1)
      pp <- pp + 
        geom_text(aes(label = sprintf('%d (%s)', ..count.., pct_format(..count.. / sum(..count..)))),
                  stat = 'count',  vjust=-0.5, size = size_counts)
    }
    if(print=="pct"){
      pct_format = scales::percent_format(accuracy = .1)
      pp <- pp + 
        geom_text(aes(label = sprintf('%s', pct_format(..count.. / sum(..count..)))),
                  stat = 'count',  vjust=-0.5, size = size_counts)
    }
    if(print=="number"){
      pp <- pp + 
        geom_text(aes(label = sprintf('%d', ..count..)),
                  stat = 'count',  vjust=-0.5, size = size_counts)
    }
  }
  else pp <- NULL
  n <- sum(!is.na(DATA.[, var]))
  if(nchar(subtitle) == 0) subtitle = paste("n:", n, " (", round(n/nrow(DATA.)*100, 1), "%)", sep="")
  pp <- pp + ylab("") + ggtitle(title, subtitle = subtitle)
  
  if(axis_angle) pp <- pp + theme(axis.text.x=element_text(angle=60,hjust=1))
  
  return(pp)
}





