asm_CosSim <- function(x,y){
     # x and y are vectors of size n where each entry is to be directly compared

     num = sum(x*y);
     
     den1 = sqrt(sum(x*x)) + 1e-10;
     
     den2 = sqrt(sum(y*y)) + 1e-10;
     
     return(num/(den1*den2))
     
}