library(latChanNet)
data(email_data)

reset_doubleChannels = function(mod){
  cur_pars = mod$get_pars()
  new_pars = cbind(cur_pars, cur_pars)
  is_pos = new_pars > 0
  n_pos = sum(is_pos)
  new_vals = runif(n_pos, max = new_pars[is_pos])
  new_pars[is_pos] = new_vals
  mod$resize(new_pars)
} 

seq_fit = function(mod, nDoubles = 3){
  mod$fit(par = T, fast_em = T)
  parList = list()
  cur_pars = mod$get_pars()
  init_cols = ncol(cur_pars)
  parList[[1]] = cur_pars
  for(i in seq_len(nDoubles) ){
    reset_doubleChannels(mod)
    mod$fit(par = T, fast_em = T)
    parList[[i + 1]] = mod$get_pars()
  }
  list_names = paste(2^(0:nDoubles), "channels")
  names(parList) = list_names
  return(parList)
}

mod1 = makeLatentModel(email_data$edgeList, 256, model = "LCN")
system.time(res1 <- mod1$fit(par = T, fast_em = T) )

mod2 = makeLatentModel(email_data$edgeList, 256, model = "LCN")
system.time(res2 <- mod2$fit(par = T, fast_em = F))