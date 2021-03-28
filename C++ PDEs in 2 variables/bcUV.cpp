result[0]=-1. + Qtt.diff("d00" ,i,j);
result[1]=-1. + Qzz.diff("d00" ,i,j);
result[2]=-1. + Qxx.diff("d00" ,i,j);
result[3]=-1. + Qyy.diff("d00" ,i,j);
result[4]=0. + Qz1.diff("d00" ,i,j);
result[5]=0. + ay.diff("d00" ,i,j);
result[6]=0. + by.diff("d00" ,i,j);
result[7]=0. + h.diff("d10" ,i,j);