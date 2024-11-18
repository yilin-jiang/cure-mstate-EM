makepre_com=function(pre_nc,pre_c,P_nc){
  pre_com=pre_c
  for(i in 1:5){
    pre_com[[i]][,-1]=pre_c[[i]][,-1]*(1-P_nc)+pre_nc[[i]][,-1]*P_nc}
  return(pre_com)
}