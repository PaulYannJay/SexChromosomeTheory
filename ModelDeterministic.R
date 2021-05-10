
### For Paper Fig 2a ##
### XY Paper Data ###

TableLinkR=data.frame(time=double(), h=double(), s=double(), u=double(), n=double(), r=double(), P=double(),
                      FXN=double(),FXI=double(),FXNm=double(),FXIm=double(),FXNf=double(),FXIf=double(),FXm=double(),FXf=double(),  
                      FYN=double(),FYI=double(), FY=double(), Wm=double(), Wf=double(),D=double(),q=double())
write.table(TableLinkR, "~/Analysis/DelSheltering/Output/TimeSimul_h_s_u1e-8_cM_PapData_2MbInv_0.80_Fl_XYsyst_Y_BC.txt", append = F, quote=F, row.names = F)
P=0.80

for (R in c(0.0, 0.001, 0.005, 0.01, 0.05, 0.5)) # Recombination between the inversion and the XY locus
{
  for (h in c(0.001, 0.01, 0.1)) #Dominance coefficient
  {
    for (s in c(0.01, 0.05, 0.1, 0.5)) #Selection coefficient
    {
      TableLinkR=data.frame(time=double(), h=double(), s=double(), u=double(), n=double(), r=double(), P=double(),
                            FXN=double(),FXI=double(),FXNm=double(),FXIm=double(),FXNf=double(),FXIf=double(),FXm=double(),FXf=double(),  
                            FYN=double(),FYI=double(), FY=double(), Wm=double(), Wf=double(),D=double(),q=double())
      u=1e-08 #Mutation rate
      n=2000000 #Inversion size
      q=((h*(1+u))/(2*(2*h -1)))*(1-sqrt(1-((4*(2*h - 1)*u)/(s*h*h*(1+u)^2)))) #mutation frequency
      r=floor(P*n*q) #Number of mutation in inversion (parameter m)
      WNI=exp(r*q*(2*h*s-s)-h*s*(n*q+r)) #Fitness of individual heterozygous for the inversion
      WNN=exp(n*q*q*(2*h*s-s)-h*s*(n*q*2)) #Fitness of individual without inversion.
      WII=exp(-s*r) #Fitness of individual homozygous for the inversion
      FXNm=1.00 #Initial frequency of non-inverted segment on the X chromosome in male
      FXIm=0.00 #Initial frequency of inversion on the X chromosome in male (=1-FXNm)
      FXNf=1.00 #Initial frequency of non-inverted segment on the X chromosome in female
      FXIf=0.00 #Initial frequency of inversion on the X chromosome in female
      FYN=0.99 #Initial frequency of non-inverted segment on the Y chromosome
      FYI=0.01 #Initial frequency of the inversion on the Y chromosome
      FY=FYN+FYI #Frequency of the Y chromosome (must be 1)
      FXm=FXNm+FXIm
      FXf=FXNf+FXIf
      FXI=(2/3)*FXIf + (1/3)*FXIm #Overall frequency of the inversion in X chromosome
      FXN=(2/3)*FXNf + (1/3)*FXNm
      Wm=FXNf*FYN*WNN + FXNf*FYI*WNI + FXIf*FYN*WNI + FXIf*FYI*WII #Male fitness
      Wf=FXNf*FXNm*WNN + FXNf*FXIm*WNI + FXNm*FXIf*WNI + FXIf*FXIm*WII #Female fitness
      D=FXI*FYN - FXN*FYI #Linkage desequilibrium
      FirstGen=c(1,h,s,u,n,R,P,FXN,FXI,FXNm,FXIm,FXNf,FXIf,FXm,FXf,FYN,FYI,FY,Wm,Wf,D,q) #State at the first generation
      TableLinkR[1,]=FirstGen 
      for (time in seq(2,10000,1)) #Simulate for 10000 generation 
      {
        FYI=(TableLinkR$FYI[time-1]*TableLinkR$FXIf[time-1]*WII +
               TableLinkR$FYI[time-1]*TableLinkR$FXNf[time-1]*WNI*(1-R) + 
               TableLinkR$FYN[time-1]*TableLinkR$FXIf[time-1]*WNI*R)/TableLinkR$Wm[time-1] #Change of FYI 
        FYN=(TableLinkR$FYN[time-1]*TableLinkR$FXNf[time-1]*WNN +
               TableLinkR$FYN[time-1]*TableLinkR$FXIf[time-1]*WNI*(1-R) + 
               TableLinkR$FYI[time-1]*TableLinkR$FXNf[time-1]*WNI*R)/TableLinkR$Wm[time-1] #Change of FYN
        
        FXIm=(TableLinkR$FXIf[time-1]*TableLinkR$FYI[time-1]*WII +
                TableLinkR$FXIf[time-1]*TableLinkR$FYN[time-1]*WNI*(1-R) + 
                TableLinkR$FXNf[time-1]*TableLinkR$FYI[time-1]*WNI*R)/TableLinkR$Wm[time-1] #Change of FXIm
        FXNm=(TableLinkR$FXNf[time-1]*TableLinkR$FYN[time-1]*WNN +
                TableLinkR$FXNf[time-1]*TableLinkR$FYI[time-1]*WNI*(1-R) + 
                TableLinkR$FXIf[time-1]*TableLinkR$FYN[time-1]*WNI*R)/TableLinkR$Wm[time-1]
        FXIf=(TableLinkR$FXIf[time-1]*TableLinkR$FXIm[time-1]*WII +
                (1/2)*(TableLinkR$FXIm[time-1]*TableLinkR$FXNf[time-1] + TableLinkR$FXIf[time-1]*TableLinkR$FXNm[time-1])*WNI)/TableLinkR$Wf[time-1]
        FXNf=(TableLinkR$FXNf[time-1]*TableLinkR$FXNm[time-1]*WNN +
                (1/2)*(TableLinkR$FXIm[time-1]*TableLinkR$FXNf[time-1] + TableLinkR$FXIf[time-1]*TableLinkR$FXNm[time-1])*WNI)/TableLinkR$Wf[time-1]
        FXI=(2/3)*FXIf + (1/3)*FXIm #To check that everything is alright
        FXN=(2/3)*FXNf + (1/3)*FXNm
        FY=FYN+FYI
        FX=FXN+FXI
        D=FXI*FYN - FXN*FYI
        Wm=FXNf*FYN*WNN + FXNf*FYI*WNI + FXIf*FYN*WNI + FXIf*FYI*WII
        Wf=FXNf*FXNm*WNN + FXNf*FXIm*WNI + FXNm*FXIf*WNI + FXIf*FXIm*WII
        Gen=c(time,h,s,u,n,R,P,FXN,FXI,FXNm,FXIm,FXNf,FXIf,FXm,FXf,FYN,FYI,FY,Wm,Wf,D,q) #State of the population at this generation
        TableLinkR[nrow(TableLinkR)+1,]=Gen
      }
      write.table(TableLinkR, "~/Analysis/DelSheltering/Output/TimeSimul_h_s_u1e-8_cM_PapData_2MbInv_0.80_Fl_XYsyst_Y_BC.txt", append = T, quote=F, row.names = F, col.names = F)
      
    }
  }
}

