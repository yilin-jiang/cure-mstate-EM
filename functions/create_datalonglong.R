create_datalonglong <- function(data_long, transition, id) {
  no.tran2 <- length(transition)
  name <- paste0("cure", transition)
  
  data_long$fromto <- paste0(data_long[,"from"], data_long[,"to"])
  
  # Initialize matrix
  transition_matrix <- matrix(0, nrow = nrow(data_long), ncol = no.tran2, 
                              dimnames = list(NULL, name))
  

  id_tran2<- which(data_long$fromto %in% transition)
  data_long2 <- cbind(data_long, transition_matrix)
  data_long2[id_tran2,]$trans <- paste0(data_long2[id_tran2,]$trans, ".2")
  
  lapply(1:no.tran2, function(i) {
    transition_idx <- which(data_long$fromto == transition[i])
    transition_matrix[transition_idx, i] <<- 1
  })

  data_long1 <- cbind(data_long, transition_matrix)
  
  data_longlong <- rbind(data_long1, data_long2)
  
  data_longlong <- unique(data_longlong)
  
  # Order by id, from, and to columns
  data_longlong <- data_longlong[order(data_longlong[, id], data_longlong$from, data_longlong$to), ]
  
  # Reset row names
  rownames(data_longlong) <- NULL
  
  return(data_longlong)
}
