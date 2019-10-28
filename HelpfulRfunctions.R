#Helpful functions, alphabetical
######################################################
#stolen from Mage
add.alpha <- function(col, alpha=1){
  if(missing(col))
    stop("Please provide a vector of colours.")
  apply(sapply(col, col2rgb)/255, 2, 
        function(x) 
          rgb(x[1], x[2], x[3], alpha=alpha))  
}
######################################################
formatter <- function(reg){
  #reg is a regression object
  
  #first determine significance
  print("Doublecheck these are the p values! If not, change indexing in function")
  print(summary(reg)$coefficients[,5])
  sigs <- vector()
  k <- 1
  for(i in 1:length(summary(reg)$coefficients[,5])){
    sigs[k] <- if(summary(reg)$coefficients[i,5] < 0.1) "*" else ""
    sigs[k] <- if(summary(reg)$coefficients[i,5] < 0.05) "**" else ""
    sigs[k] <- if(summary(reg)$coefficients[i,5] < 0.01) "***" else ""
    k <- k + 1
  }
  
  paste(round(summary(reg)$coefficients[,1],1), 
        sigs,
        " (", 
        round(summary(reg)$coefficients[,1] - 1.96*summary(reg)$coefficients[,2],1),
        ",",
        round(summary(reg)$coefficients[,1] + 1.96*summary(reg)$coefficients[,2],1),
        ")",
        sep=""
  )
}
######################################################
plotR <- function(y, x1, x2, ptcol){
  add.alpha <- function(col, alpha=1){
    if(missing(col))
      stop("Please provide a vector of colours.")
    apply(sapply(col, col2rgb)/255, 2, 
          function(x) 
            rgb(x[1], x[2], x[3], alpha=alpha))  
  }
  stripchart(y~x1 + x2,
             vertical = T,
             xaxt="n",
             yaxt="n",
             pch = 21,
             cex = 1.7,
             ylab = "",
             xlab = "",
             col =  add.alpha(ptcol, 0.8),
             bg = add.alpha(ptcol, 0.4),
             method = "jitter",
             frame.plot = F)
  boxplot(y~x1 + x2,
          add = T,
          col = add.alpha(ptcol, 0.4),
          frame.plot = F,
          yaxt = "n",
          xaxt = "n",
          ylab = "",
          xlab = "",
          outline = F
  )
}

######################################################