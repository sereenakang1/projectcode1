# Import pandas as pd
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt 


##################################################################################################
#Loading Pu Seq data

#chromosome 1 
fddata = np.load("chr1_f_d.npy")#plot with forward epsilon
fedata = np.load("chr1_f_e.npy")
ordata = np.load("chr1_orieff.npy")

redata = np.load("chr1_r_e.npy")
rddata = np.load("chr1_r_d.npy")



#chromosome 2 
fddata2 = np.load("chr2_f_d.npy")#plot with forward epsilon
fedata2 = np.load("chr2_f_e.npy")
ordata2 = np.load("chr2_orieff.npy")

redata2 = np.load("chr2_r_e.npy")
rddata2 = np.load("chr2_r_d.npy")


#chromosome 3 


fddata3 = np.load("chr3_f_d.npy")#plot with forward epsilon
fedata3 = np.load("chr3_f_e.npy")
ordata3 = np.load("chr3_orieff.npy")

redata3 = np.load("chr3_r_e.npy")
rddata3 = np.load("chr3_r_d.npy")

#Creating x axes for Fe Data

x1 = np.arange(150, len(fedata)*300,300) 

x2 = np.arange(150, len(fedata2)*300,300) 

x3 = np.arange(150, len(fedata3)*300,300) 




#########################################################################

#I write ef for epsilon forward, dr for delta reverse.
average = 0.5*(fedata+rddata) #average of epsilon forwad and delta reverse

d_rt = np.gradient(average) # this computes the "gradient" (the slope)
cl_drt = np.clip(d_rt,0,None) #This sets to 0 all the places the had a negative gradient (downward slope).



average2 = 0.5*(fedata2+rddata2) #average of epsilon forwad and delta reverse

d_rt2 = np.gradient(average2) # this computes the "gradient" (the slope)
cl_drt2 = np.clip(d_rt2,0,None)



average3 = 0.5*(fedata3+rddata3) #average of epsilon forwad and delta reverse

d_rt3 = np.gradient(average3) # this computes the "gradient" (the slope)
cl_drt3 = np.clip(d_rt3,0,None)


plt.plot(x2,cl_drt2)
plt.xlim(0,20000)
plt.ylabel("Origin Efficiency")
plt.xlabel("Chromosome II position")
plt.title("Average fraction of right forks - origin efficiency")
plt.show()


plt.plot(x3,cl_drt3)
plt.xlim(0,20000)
plt.ylabel("Origin Efficiency")
plt.xlabel("Chromosome III position")
plt.title("Average fraction of right forks - origin efficiency")
plt.show()


#Plotting replication direction
#when writing up: explain how replication direction was calculated

direction = (2*fedata)-1

plt.plot(x1,direction)
plt.xlim(0, 20000)
plt.xlabel("Chromosome I position")
plt.ylabel("Replication Direction")
plt.title("Replication Direction")
plt.show()


direction1 = (2*fedata2)-1

plt.plot(x2,direction1)
plt.xlim(0, 20000)
plt.xlabel("Chromosome II position")
plt.ylabel("Replication Direction")
plt.title("Replication Direction")
plt.show()





direction2 = (2*fedata3)-1

plt.plot(x3,direction2)
plt.xlim(0, 20000)
plt.xlabel("Chromosome III position")
plt.ylabel("Replication Direction")
plt.title("Replication Direction")
plt.show()


#Plotting Puseq data


#Chromosome 1





fig, (ax1,ax2) = plt.subplots(nrows=2, ncols=1, sharex = True)

ax1.plot(x1, fedata, label='Delta Polymerase')
ax1.plot(x1,fddata, label = "Epsilon Polymerase")
ax2.plot(x1,redata, x1, rddata)#epsilion = orange
plt.xlim(0, 20000)
ax1.set_title("Forward")
ax2.set_title("Reverse")
plt.xlabel("Chromosome I position (basepairs)")
fig.text(0.06, 0.5, 'Ratio', ha='center', va='center', rotation='vertical')
plt.suptitle("Polymerase Usage")

ax1.legend(bbox_to_anchor=(1.05, 0.1), loc='upper left', borderaxespad=0.)
plt.show()


#Chromosome 2

fig, (ax3,ax4) = plt.subplots(nrows=2, ncols=1, sharex = True)

ax3.plot(x2,fedata2, label='Delta Polymerase' )
ax3.plot(x2,fddata2, label='Epsilon Polymerase')
ax3.legend(bbox_to_anchor=(1.05, 0.1), loc='upper left', borderaxespad=0.)
ax4.plot(x2,redata2, x2, rddata2)
plt.xlim(0, 20000)
ax3.set_title("Forward")
ax4.set_title("Reverse")
plt.xlabel("Chromosome II position (kb)")
fig.text(0.06, 0.5, 'Ratio', ha='center', va='center', rotation='vertical')
plt.suptitle("Polymerase Usage")
plt.show()


#Chromosome 3 

fig, (ax3,ax4) = plt.subplots(nrows=2, ncols=1, sharex = True)

ax3.plot(x3, fedata3, label='Delta Polymerase')
ax3.plot(x3, fddata3, label='Epsilon Polymerase')
ax4.plot(x3,redata3,x3, rddata3)
ax3.legend(bbox_to_anchor=(1.05, 0.1), loc='upper left', borderaxespad=0.)
plt.xlim(0, 20000)
ax3.set_title("Forward")
ax4.set_title("Reverse")
plt.xlabel("Chromosome III position (kb)")
fig.text(0.06, 0.5, 'Ratio', ha='center', va='center', rotation='vertical')
plt.suptitle("Polymerase Usage")
plt.show()

fig, (ax7,ax8,ax9) = plt.subplots(nrows =3, ncols =1, sharex = True)
ax7.plot(x1, ordata)
ax8.plot(x2, ordata2)
ax9.plot(x3, ordata3)
plt.xlabel("Chromosome position")
fig.text(0.06, 0.5, 'Origin Efficiency', ha='center', va='center', rotation='vertical')
plt.xlim(0,20000)
plt.show()

plt.plot(x1,cl_drt)
plt.xlim(0,20000)
plt.ylabel("Origin Efficiency")
plt.xlabel("Chromosome I position")
plt.title("Average fraction of right forks - origin efficiency")
############################################################################
#Loading transcription table + selecting specific columns

transcript = pd.read_csv('pombe_transcription.csv')

df = pd.DataFrame(transcript) 


chr1 = df.loc[df['chromosomea']==1] #gets data for each specfic chromosome 
chr2 = df.loc[df['chromosomea']==2]
chr3 = df.loc[df['chromosomea']==3]

rna1 = chr1[['MM.mRNA.cpc']] # Get RNA levels 


#calculating base sizes
Object = pd.DataFrame(chr1, columns = ('ORF.start','ORF.end' ))
differenceFrame = Object.diff(axis=1) 
differenceFrame.rename(columns={'ORF.end':'Base Size'}, inplace=True) # contains column with difference between ORF start and end. 

differenceFrame.drop('ORF.start', axis=1, inplace=True)
########################################################################################################
#Splitting by direction 

strand = [chr1['coding.strand']]
plusdirection = chr1.loc[chr1['coding.strand']== '+']
minusdirection = chr1.loc[chr1['coding.strand']== '-']

startbinplus = [ x // 300 for x in plusdirection["ORF.start"]  ]
endbinplus = [ x // 300 for x in plusdirection["ORF.end"] ]

Z = np.linspace(0,1,100) 
resultplus = [] 
for s,e in zip(startbinplus,endbinplus):
    Z1 = np.linspace(0,1,(e-s)+1) 
    resultplus.append(np.interp(Z, Z1, fedata[s:e+1]))

resnplus = np.array(resultplus) 
meanplus = np.nanmean(resnplus,axis = 0)




startbinminus = [ x // 300 for x in minusdirection["ORF.start"]  ]
endbinminus = [ x // 300 for x in minusdirection["ORF.end"] ]

Z = np.linspace(0,1,100) 
resultminus = [] 
for s,e in zip(startbinminus,endbinminus):
    Z1 = np.linspace(0,1,(e-s)+1) 
    resultminus.append(np.interp(Z, Z1, redata[s:e+1]))

resnpminus = np.array(resultminus) 
meanminus = np.nanmean(resnpminus,axis = 0)

fig, (ax16,ax17) = plt.subplots(2, sharex = True, sharey = True)
ax16.plot(meanplus)
ax17.plot(meanminus)

fig.text(0.02, 0.5, 'Average Forward Epsilon usage', ha='center', va='center', rotation='vertical')
plt.xlabel("Gene bins")
plt.suptitle("Genes split by direction (chromosome 1) ")

#########################################################################################






############################################################################################


plusdirection3 = chr3.loc[chr3['coding.strand']== '+']
minusdirection3 = chr3.loc[chr3['coding.strand']== '-']


##########################################################################################################
#Sorting genes into quartiles of transcription activity and plotting


q25 = df[ df["MM.mRNA.cpc"] < df["MM.mRNA.cpc"].quantile(.25) ] #From percentile 0% to 25%
q25chr1 = q25[q25['chromosomea']==1]

startbin25 = [ x // 300 for x in q25chr1["ORF.start"]  ]
endbin25 = [ x // 300 for x in q25chr1["ORF.end"] ]

A = np.linspace(0,1,100) 
result25 = [] 
for s,e in zip(startbin25,endbin25):
    
    A1 = np.linspace(0,1,len(fedata[s:e+1])) 
    result25.append(np.interp(A, A1, fedata[s:e+1]))

resnp25 = np.array(result25) 
mean25 = np.nanmean(resnp25,axis = 0)








q50 = df[ (df["MM.mRNA.cpc"].quantile(.25) < df["MM.mRNA.cpc"]) & (df["MM.mRNA.cpc"] < df["MM.mRNA.cpc"].quantile(.5))] #25% to 50%
q50chr1 = q50[q50['chromosomea']==1]

startbin50 = [ x // 300 for x in q50chr1["ORF.start"]  ]
endbin50 = [ x // 300 for x in q50chr1["ORF.end"] ]

B = np.linspace(0,1,100) 
result50 = [] 
for s,e in zip(startbin50,endbin50):
    
    B1 = np.linspace(0,1,len(fedata[s:e+1])) 
    result50.append(np.interp(B, B1, fedata[s:e+1]))

resnp50 = np.array(result50) 
mean50 = np.nanmean(resnp50,axis = 0)




q75 = df[ (df["MM.mRNA.cpc"].quantile(.5) < df["MM.mRNA.cpc"]) & (df["MM.mRNA.cpc"] < df["MM.mRNA.cpc"].quantile(.75))] #50% to 75%
q75chr1 = q75[q75['chromosomea']==1]

startbin75 = [ x // 300 for x in q75chr1["ORF.start"]  ]
endbin75 = [ x // 300 for x in q75chr1["ORF.end"] ]

C = np.linspace(0,1,100) 
result75 = [] 
for s,e in zip(startbin75,endbin75):
    
    C1 = np.linspace(0,1,len(fedata[s:e+1])) 
    result75.append(np.interp(C, C1, fedata[s:e+1]))

resnp75 = np.array(result75) 
mean75 = np.nanmean(resnp75,axis = 0)





q100 = df[ df["MM.mRNA.cpc"].quantile(.75) < df["MM.mRNA.cpc"] ]

q100chr1 = q100[q100['chromosomea']==1]

startbin100 = [ x // 300 for x in q100chr1["ORF.start"]  ]
endbin100 = [ x // 300 for x in q100chr1["ORF.end"] ]

C = np.linspace(0,1,100) 
result100 = [] 
for s,e in zip(startbin100,endbin100):
    
    C1 = np.linspace(0,1,len(fedata[s:e+1])) 
    result100.append(np.interp(C, C1, fedata[s:e+1]))

resnp100 = np.array(result100) 
mean100 = np.nanmean(resnp100,axis = 0)



#plotting transcription levels


fig, axs = plt.subplots(2, 2, sharex = True, sharey = True)
axs[0, 0].plot(mean25)
fig.text(0.02, 0.5, 'Average Forward Epsilon usage', ha='center', va='center', rotation='vertical')
plt.xlabel('Gene Bins - Chromosome 1', horizontalalignment='center',position=(0.0001,30))
plt.suptitle("Genes split into quantiles from lowest to highest RNA Levels")
axs[0, 0].set_title('0-25%')
axs[0, 1].plot(mean50, 'tab:orange')
axs[0, 1].set_title('25%-50%')
axs[1, 0].plot(mean75, 'tab:green')
axs[1, 0].set_title('50%-75%')
axs[1, 1].plot(mean100, 'tab:red')
axs[1, 1].set_title('75%-100%')
plt.show()
#########################################################################################
q25 = df[ df["MM.mRNA.cpc"] < df["MM.mRNA.cpc"].quantile(.25) ] #From percentile 0% to 25%
q25chr2 = q25[q25['chromosomea']==2]

startbin25 = [ x // 300 for x in q25chr2["ORF.start"]  ]
endbin25 = [ x // 300 for x in q25chr2["ORF.end"] ]

A = np.linspace(0,1,100) 
result25 = [] 
for s,e in zip(startbin25,endbin25):
    
    A1 = np.linspace(0,1,len(fedata2[s:e+1])) 
    result25.append(np.interp(A, A1, fedata2[s:e+1]))

resnp25 = np.array(result25) 
mean25 = np.nanmean(resnp25,axis = 0)








q50 = df[ (df["MM.mRNA.cpc"].quantile(.25) < df["MM.mRNA.cpc"]) & (df["MM.mRNA.cpc"] < df["MM.mRNA.cpc"].quantile(.5))] #25% to 50%
q50chr2 = q50[q50['chromosomea']==2]

startbin50 = [ x // 300 for x in q50chr2["ORF.start"]  ]
endbin50 = [ x // 300 for x in q50chr2["ORF.end"] ]

B = np.linspace(0,1,100) 
result50 = [] 
for s,e in zip(startbin50,endbin50):
    
    B1 = np.linspace(0,1,len(fedata2[s:e+1])) 
    result50.append(np.interp(B, B1, fedata2[s:e+1]))

resnp50 = np.array(result50) 
mean50 = np.nanmean(resnp50,axis = 0)




q75 = df[ (df["MM.mRNA.cpc"].quantile(.5) < df["MM.mRNA.cpc"]) & (df["MM.mRNA.cpc"] < df["MM.mRNA.cpc"].quantile(.75))] #50% to 75%
q75chr2 = q75[q75['chromosomea']==2]

startbin75 = [ x // 300 for x in q75chr2["ORF.start"]  ]
endbin75 = [ x // 300 for x in q75chr2["ORF.end"] ]

C = np.linspace(0,1,100) 
result75 = [] 
for s,e in zip(startbin75,endbin75):
    
    C1 = np.linspace(0,1,len(fedata2[s:e+1])) 
    result75.append(np.interp(C, C1, fedata2[s:e+1]))

resnp75 = np.array(result75) 
mean75 = np.nanmean(resnp75,axis = 0)





q100 = df[ df["MM.mRNA.cpc"].quantile(.75) < df["MM.mRNA.cpc"] ]

q100chr2 = q100[q100['chromosomea']==2]

startbin100 = [ x // 300 for x in q100chr2["ORF.start"]  ]
endbin100 = [ x // 300 for x in q100chr2["ORF.end"] ]

C = np.linspace(0,1,100) 
result100 = [] 
for s,e in zip(startbin100,endbin100):
    
    C1 = np.linspace(0,1,len(fedata2[s:e+1])) 
    result100.append(np.interp(C, C1, fedata2[s:e+1]))

resnp100 = np.array(result100) 
mean100 = np.nanmean(resnp100,axis = 0)



#plotting transcription levels


fig, axs = plt.subplots(2, 2, sharex = True, sharey = True)
axs[0, 0].plot(mean25)
fig.text(0.02, 0.5, 'Average Forward Epsilon usage', ha='center', va='center', rotation='vertical')
plt.xlabel('Gene Bins - Chromosome 2', horizontalalignment='center',position=(0.0001,30))
plt.suptitle("Genes split into quantiles from lowest to highest RNA Levels")
axs[0, 0].set_title('0-25%')
axs[0, 1].plot(mean50, 'tab:orange')
axs[0, 1].set_title('25%-50%')
axs[1, 0].plot(mean75, 'tab:green')
axs[1, 0].set_title('50%-75%')
axs[1, 1].plot(mean100, 'tab:red')
axs[1, 1].set_title('75%-100%')
plt.show()

################################################################################################

q25 = df[ df["MM.mRNA.cpc"] < df["MM.mRNA.cpc"].quantile(.25) ] #From percentile 0% to 25%
q25chr3 = q25[q25['chromosomea']==3]

startbin25 = [ x // 300 for x in q25chr3["ORF.start"]  ]
endbin25 = [ x // 300 for x in q25chr3["ORF.end"] ]

A = np.linspace(0,1,100) 
result25 = [] 
for s,e in zip(startbin25,endbin25):
    
    A1 = np.linspace(0,1,len(fedata3[s:e+1])) 
    result25.append(np.interp(A, A1, fedata3[s:e+1]))

resnp25 = np.array(result25) 
mean25 = np.nanmean(resnp25,axis = 0)








q50 = df[ (df["MM.mRNA.cpc"].quantile(.25) < df["MM.mRNA.cpc"]) & (df["MM.mRNA.cpc"] < df["MM.mRNA.cpc"].quantile(.5))] #25% to 50%
q50chr3 = q50[q50['chromosomea']==3]

startbin50 = [ x // 300 for x in q50chr3["ORF.start"]  ]
endbin50 = [ x // 300 for x in q50chr3["ORF.end"] ]

B = np.linspace(0,1,100) 
result50 = [] 
for s,e in zip(startbin50,endbin50):
    
    B1 = np.linspace(0,1,len(fedata3[s:e+1])) 
    result50.append(np.interp(B, B1, fedata3[s:e+1]))

resnp50 = np.array(result50) 
mean50 = np.nanmean(resnp50,axis = 0)




q75 = df[ (df["MM.mRNA.cpc"].quantile(.5) < df["MM.mRNA.cpc"]) & (df["MM.mRNA.cpc"] < df["MM.mRNA.cpc"].quantile(.75))] #50% to 75%
q75chr3 = q75[q75['chromosomea']==3]

startbin75 = [ x // 300 for x in q75chr3["ORF.start"]  ]
endbin75 = [ x // 300 for x in q75chr3["ORF.end"] ]

C = np.linspace(0,1,100) 
result75 = [] 
for s,e in zip(startbin75,endbin75):
    
    C1 = np.linspace(0,1,len(fedata3[s:e+1])) 
    result75.append(np.interp(C, C1, fedata3[s:e+1]))

resnp75 = np.array(result75) 
mean75 = np.nanmean(resnp75,axis = 0)





q100 = df[ df["MM.mRNA.cpc"].quantile(.75) < df["MM.mRNA.cpc"] ]

q100chr3 = q100[q100['chromosomea']==3]

startbin100 = [ x // 300 for x in q100chr3["ORF.start"]  ]
endbin100 = [ x // 300 for x in q100chr3["ORF.end"] ]

C = np.linspace(0,1,100) 
result100 = [] 
for s,e in zip(startbin100,endbin100):
    
    C1 = np.linspace(0,1,len(fedata3[s:e+1])) 
    result100.append(np.interp(C, C1, fedata3[s:e+1]))

resnp100 = np.array(result100) 
mean100 = np.nanmean(resnp100,axis = 0)



#plotting transcription levels


fig, axs = plt.subplots(2, 2, sharex = True, sharey = True)
axs[0, 0].plot(mean25)
fig.text(0.02, 0.5, 'Average Forward Epsilon usage', ha='center', va='center', rotation='vertical')
plt.xlabel('Gene Bins - Chromosome 3', horizontalalignment='center',position=(0.0001,30))
plt.suptitle("Genes split into quantiles from lowest to highest RNA Levels")
axs[0, 0].set_title('0-25%')
axs[0, 1].plot(mean50, 'tab:orange')
axs[0, 1].set_title('25%-50%')
axs[1, 0].plot(mean75, 'tab:green')
axs[1, 0].set_title('50%-75%')
axs[1, 1].plot(mean100, 'tab:red')
axs[1, 1].set_title('75%-100%')
plt.show()
















#######################################################################################################

#quantile plot - 1.5kb before and after, split between forward/reverse strands

#forward strand - before and after




plus25 = q25chr1.loc[q25chr1['coding.strand']== '+']
plus50 = q50chr1.loc[q50chr1['coding.strand']== '+']
plus75 = q75chr1.loc[q75chr1['coding.strand']== '+']
plus100 = q100chr1.loc[q100chr1['coding.strand']== '+']




startbin25plus = [ x // 300 for x in plus25["ORF.start"]  ]
endbin25plus = [ x // 300 for x in plus25["ORF.end"] ]

A = np.linspace(0,1,100) 
result25plus = [] 
resultbefore25 = []
resultafter25 = []
for s,e in zip(startbin25plus,endbin25plus):
    
    A1 = np.linspace(0,1,len(fedata[s:e+1])) 
    result25plus.append(np.interp(A, A1, fedata[s:e+1]))
    resultbefore25.append(fedata[s-5:s])
    resultafter25.append(fedata[e:e+5])
resnp25plus = np.array(result25plus) 
mean25plus = np.nanmean(resnp25,axis = 0)


startbin50plus = [ x // 300 for x in plus50["ORF.start"]  ]
endbin50plus = [ x // 300 for x in plus50["ORF.end"] ]

A = np.linspace(0,1,100) 
result50plus = [] 
resultbefore50 = []
resultafter50 = []
for s,e in zip(startbin50plus,endbin50plus):
    
    A1 = np.linspace(0,1,len(fedata[s:e+1])) 
    result50plus.append(np.interp(A, A1, fedata[s:e+1]))
    resultbefore50.append(fedata[s-5:s])
    resultafter50.append(fedata[e:e+5])
resnp50plus = np.array(result50plus) 
mean50plus = np.nanmean(resnp50plus,axis = 0)



startbin75plus = [ x // 300 for x in plus75["ORF.start"]  ]
endbin75plus = [ x // 300 for x in plus75["ORF.end"] ]

A = np.linspace(0,1,100) 
result75plus = [] 
resultbefore75 = []
resultafter75 = []
for s,e in zip(startbin75plus,endbin75plus):
    
    A1 = np.linspace(0,1,len(fedata[s:e+1])) 
    result75plus.append(np.interp(A, A1, fedata[s:e+1]))
    resultbefore75.append(fedata[s-5:s])
    resultafter75.append(fedata[e:e+5])
resnp75plus = np.array(result75plus) 
mean75plus = np.nanmean(resnp75plus,axis = 0)



startbin100plus = [ x // 300 for x in plus100["ORF.start"]  ]
endbin100plus = [ x // 300 for x in plus100["ORF.end"] ]

A = np.linspace(0,1,100) 
result100plus = [] 
resultbefore100 = []
resultafter100 = []
for s,e in zip(startbin100plus,endbin100plus):
    
    A1 = np.linspace(0,1,len(fedata[s:e+1])) 
    result100plus.append(np.interp(A, A1, fedata[s:e+1]))
    resultbefore100.append(fedata[s-5:s])
    resultafter100.append(fedata[e:e+5])
resnp100plus = np.array(result100plus) 
mean100plus = np.nanmean(resnp100plus,axis = 0)



fig, axs = plt.subplots(2, 2, sharex = True, sharey = True)
axs[0, 0].plot(mean25plus)
fig.text(0.02, 0.5, 'Average FE usage, + strand', ha='center', va='center', rotation='vertical')
plt.xlabel("Gene bins - Chromosome 1", horizontalalignment='center', position = (0.01, 25))
plt.suptitle("Genes split into quantiles from lowest to highest RNA Levels")
axs[0, 0].set_title('0-25%')
axs[0, 1].plot(mean50plus, 'tab:orange')
axs[0, 1].set_title('25%-50%')
axs[1, 0].plot(mean75plus, 'tab:green')
axs[1, 0].set_title('50%-75%')
axs[1, 1].plot(mean100plus, 'tab:red')
axs[1, 1].set_title('75%-100%')
plt.show()


#minus strand quantiles

minus25 = q25chr1.loc[q25chr1['coding.strand']== '-']
minus50 = q50chr1.loc[q50chr1['coding.strand']== '-']
minus75 = q75chr1.loc[q75chr1['coding.strand']== '-']
minus100 = q100chr1.loc[q100chr1['coding.strand']== '-']




startbin25minus = [ x // 300 for x in minus25["ORF.start"]  ]
endbin25minus = [ x // 300 for x in minus25["ORF.end"] ]

A = np.linspace(0,1,100) 
result25 = [] 
resultbefore25 = []
resultafter25 = []
for s,e in zip(startbin25minus,endbin25minus):
    
    A1 = np.linspace(0,1,len(redata[s:e+1])) 
    result25.append(np.interp(A, A1, redata[s:e+1][::-1]))
    resultbefore25.append(redata[s-5:s][::-1])
    resultafter25.append(redata[e:e+5][::-1])
resnp25 = np.array(result25) 
mean25 = np.nanmean(resnp25,axis = 0)


startbin50minus = [ x // 300 for x in minus50["ORF.start"]  ]
endbin50minus = [ x // 300 for x in minus50["ORF.end"] ]

A = np.linspace(0,1,100) 
result50 = [] 
resultbefore50 = []
resultafter50 = []
for s,e in zip(startbin50minus,endbin50minus):
    
    A1 = np.linspace(0,1,len(redata[s:e+1])) 
    result50.append(np.interp(A, A1, redata[s:e+1][::-1]))
    resultbefore50.append(redata[s-5:s][::-1])
    resultafter50.append(redata[e:e+5][::-1])
resnp50 = np.array(result50) 
mean50 = np.nanmean(resnp50,axis = 0)



startbin75minus = [ x // 300 for x in minus75["ORF.start"]  ]
endbin75minus = [ x // 300 for x in minus75["ORF.end"] ]

A = np.linspace(0,1,100) 
result75 = [] 
resultbefore75 = []
resultafter75 = []
for s,e in zip(startbin75minus,endbin75minus):
    
    A1 = np.linspace(0,1,len(redata[s:e+1])) 
    result75.append(np.interp(A, A1, redata[s:e+1][::-1]))
    resultbefore75.append(redata[s-5:s][::-1])
    resultafter75.append(redata[e:e+5][::-1])
resnp75 = np.array(result75) 
mean75 = np.nanmean(resnp75,axis = 0)



startbin100minus = [ x // 300 for x in minus100["ORF.start"]  ]
endbin100minus = [ x // 300 for x in minus100["ORF.end"] ]

A = np.linspace(0,1,100) 
result100 = [] 
resultbefore100 = []
resultafter100 = []
for s,e in zip(startbin100minus,endbin100minus):
    
    A1 = np.linspace(0,1,len(redata[s:e+1])) 
    result100.append(np.interp(A, A1, redata[s:e+1][::-1]))
    resultbefore100.append(redata[s-5:s][::-1])
    resultafter100.append(redata[e:e+5][::-1])
resnp100 = np.array(result100) 
mean100 = np.nanmean(resnp100,axis = 0)


#reverse epsilon usage 

fig, axs = plt.subplots(2, 2, sharex = True, sharey = True)
axs[0, 0].plot(mean25)
fig.text(0.02, 0.5, 'Average RE usage, - strand', ha='center', va='center', rotation='vertical')
plt.xlabel("Gene bins - Chromosome 1", horizontalalignment='center', position = (0.01, 25))
plt.suptitle("Genes split into quantiles from lowest to highest RNA Levels")
axs[0, 0].set_title('0-25%')
axs[0, 1].plot(mean50, 'tab:orange')
axs[0, 1].set_title('25%-50%')
axs[1, 0].plot(mean75, 'tab:green')
axs[1, 0].set_title('50%-75%')
axs[1, 1].plot(mean100, 'tab:red')
axs[1, 1].set_title('75%-100%')
plt.show()
###################################################################################################
plus25 = q25chr2.loc[q25chr2['coding.strand']== '+']
plus50 = q50chr2.loc[q50chr2['coding.strand']== '+']
plus75 = q75chr2.loc[q75chr2['coding.strand']== '+']
plus100 = q100chr2.loc[q100chr2['coding.strand']== '+']




startbin25plus = [ x // 300 for x in plus25["ORF.start"]  ]
endbin25plus = [ x // 300 for x in plus25["ORF.end"] ]

A = np.linspace(0,1,100) 
result25plus = [] 
resultbefore25 = []
resultafter25 = []
for s,e in zip(startbin25plus,endbin25plus):
    
    A1 = np.linspace(0,1,len(fedata2[s:e+1])) 
    result25plus.append(np.interp(A, A1, fedata2[s:e+1]))
    resultbefore25.append(fedata2[s-5:s])
    resultafter25.append(fedata2[e:e+5])
resnp25plus = np.array(result25plus) 
mean25plus = np.nanmean(resnp25,axis = 0)


startbin50plus = [ x // 300 for x in plus50["ORF.start"]  ]
endbin50plus = [ x // 300 for x in plus50["ORF.end"] ]

A = np.linspace(0,1,100) 
result50plus = [] 
resultbefore50 = []
resultafter50 = []
for s,e in zip(startbin50plus,endbin50plus):
    
    A1 = np.linspace(0,1,len(fedata2[s:e+1])) 
    result50plus.append(np.interp(A, A1, fedata2[s:e+1]))
    resultbefore50.append(fedata2[s-5:s])
    resultafter50.append(fedata2[e:e+5])
resnp50plus = np.array(result50plus) 
mean50plus = np.nanmean(resnp50plus,axis = 0)



startbin75plus = [ x // 300 for x in plus75["ORF.start"]  ]
endbin75plus = [ x // 300 for x in plus75["ORF.end"] ]

A = np.linspace(0,1,100) 
result75plus = [] 
resultbefore75 = []
resultafter75 = []
for s,e in zip(startbin75plus,endbin75plus):
    
    A1 = np.linspace(0,1,len(fedata2[s:e+1])) 
    result75plus.append(np.interp(A, A1, fedata2[s:e+1]))
    resultbefore75.append(fedata2[s-5:s])
    resultafter75.append(fedata2[e:e+5])
resnp75plus = np.array(result75plus) 
mean75plus = np.nanmean(resnp75plus,axis = 0)



startbin100plus = [ x // 300 for x in plus100["ORF.start"]  ]
endbin100plus = [ x // 300 for x in plus100["ORF.end"] ]

A = np.linspace(0,1,100) 
result100plus = [] 
resultbefore100 = []
resultafter100 = []
for s,e in zip(startbin100plus,endbin100plus):
    
    A1 = np.linspace(0,1,len(fedata2[s:e+1])) 
    result100plus.append(np.interp(A, A1, fedata2[s:e+1]))
    resultbefore100.append(fedata2[s-5:s])
    resultafter100.append(fedata2[e:e+5])
resnp100plus = np.array(result100plus) 
mean100plus = np.nanmean(resnp100plus,axis = 0)



fig, axs = plt.subplots(2, 2, sharex = True, sharey = True)
axs[0, 0].plot(mean25plus)
fig.text(0.02, 0.5, 'Average FE usage, + strand', ha='center', va='center', rotation='vertical')
plt.xlabel("Gene bins - Chromosome 2", horizontalalignment='center', position = (0.01, 25))
plt.suptitle("Genes split into quantiles from lowest to highest RNA Levels")
axs[0, 0].set_title('0-25%')
axs[0, 1].plot(mean50plus, 'tab:orange')
axs[0, 1].set_title('25%-50%')
axs[1, 0].plot(mean75plus, 'tab:green')
axs[1, 0].set_title('50%-75%')
axs[1, 1].plot(mean100plus, 'tab:red')
axs[1, 1].set_title('75%-100%')
plt.show()


#minus strand quantiles

minus25 = q25chr2.loc[q25chr2['coding.strand']== '-']
minus50 = q50chr2.loc[q50chr2['coding.strand']== '-']
minus75 = q75chr2.loc[q75chr2['coding.strand']== '-']
minus100 = q100chr2.loc[q100chr2['coding.strand']== '-']




startbin25minus = [ x // 300 for x in minus25["ORF.start"]  ]
endbin25minus = [ x // 300 for x in minus25["ORF.end"] ]

A = np.linspace(0,1,100) 
result25 = [] 
resultbefore25 = []
resultafter25 = []
for s,e in zip(startbin25minus,endbin25minus):
    
    A1 = np.linspace(0,1,len(redata2[s:e+1])) 
    result25.append(np.interp(A, A1, redata2[s:e+1][::-1]))
    resultbefore25.append(redata2[s-5:s][::-1])
    resultafter25.append(redata2[e:e+5][::-1])
resnp25 = np.array(result25) 
mean25 = np.nanmean(resnp25,axis = 0)


startbin50minus = [ x // 300 for x in minus50["ORF.start"]  ]
endbin50minus = [ x // 300 for x in minus50["ORF.end"] ]

A = np.linspace(0,1,100) 
result50 = [] 
resultbefore50 = []
resultafter50 = []
for s,e in zip(startbin50minus,endbin50minus):
    
    A1 = np.linspace(0,1,len(redata2[s:e+1])) 
    result50.append(np.interp(A, A1, redata2[s:e+1][::-1]))
    resultbefore50.append(redata2[s-5:s][::-1])
    resultafter50.append(redata2[e:e+5][::-1])
resnp50 = np.array(result50) 
mean50 = np.nanmean(resnp50,axis = 0)



startbin75minus = [ x // 300 for x in minus75["ORF.start"]  ]
endbin75minus = [ x // 300 for x in minus75["ORF.end"] ]

A = np.linspace(0,1,100) 
result75 = [] 
resultbefore75 = []
resultafter75 = []
for s,e in zip(startbin75minus,endbin75minus):
    
    A1 = np.linspace(0,1,len(redata2[s:e+1])) 
    result75.append(np.interp(A, A1, redata2[s:e+1][::-1]))
    resultbefore75.append(redata2[s-5:s][::-1])
    resultafter75.append(redata2[e:e+5][::-1])
resnp75 = np.array(result75) 
mean75 = np.nanmean(resnp75,axis = 0)



startbin100minus = [ x // 300 for x in minus100["ORF.start"]  ]
endbin100minus = [ x // 300 for x in minus100["ORF.end"] ]

A = np.linspace(0,1,100) 
result100 = [] 
resultbefore100 = []
resultafter100 = []
for s,e in zip(startbin100minus,endbin100minus):
    
    A1 = np.linspace(0,1,len(redata2[s:e+1])) 
    result100.append(np.interp(A, A1, redata2[s:e+1][::-1]))
    resultbefore100.append(redata2[s-5:s][::-1])
    resultafter100.append(redata2[e:e+5][::-1])
resnp100 = np.array(result100) 
mean100 = np.nanmean(resnp100,axis = 0)


#reverse epsilon usage 

fig, axs = plt.subplots(2, 2, sharex = True, sharey = True)
axs[0, 0].plot(mean25)
fig.text(0.02, 0.5, 'Average RE usage, - strand', ha='center', va='center', rotation='vertical')
plt.xlabel("Gene bins - Chromosome 2", horizontalalignment='center', position = (0.01, 25))
plt.suptitle("Genes split into quantiles from lowest to highest RNA Levels")
axs[0, 0].set_title('0-25%')
axs[0, 1].plot(mean50, 'tab:orange')
axs[0, 1].set_title('25%-50%')
axs[1, 0].plot(mean75, 'tab:green')
axs[1, 0].set_title('50%-75%')
axs[1, 1].plot(mean100, 'tab:red')
axs[1, 1].set_title('75%-100%')
plt.show()
#########################################################################################
plus25 = q25chr3.loc[q25chr3['coding.strand']== '+']
plus50 = q50chr3.loc[q50chr3['coding.strand']== '+']
plus75 = q75chr3.loc[q75chr3['coding.strand']== '+']
plus100 = q100chr3.loc[q100chr3['coding.strand']== '+']




startbin25plus = [ x // 300 for x in plus25["ORF.start"]  ]
endbin25plus = [ x // 300 for x in plus25["ORF.end"] ]

A = np.linspace(0,1,100) 
result25plus = [] 
resultbefore25 = []
resultafter25 = []
for s,e in zip(startbin25plus,endbin25plus):
    
    A1 = np.linspace(0,1,len(fedata3[s:e+1])) 
    result25plus.append(np.interp(A, A1, fedata3[s:e+1]))
    resultbefore25.append(fedata3[s-5:s])
    resultafter25.append(fedata3[e:e+5])
resnp25plus = np.array(result25plus) 
mean25plus = np.nanmean(resnp25,axis = 0)


startbin50plus = [ x // 300 for x in plus50["ORF.start"]  ]
endbin50plus = [ x // 300 for x in plus50["ORF.end"] ]

A = np.linspace(0,1,100) 
result50plus = [] 
resultbefore50 = []
resultafter50 = []
for s,e in zip(startbin50plus,endbin50plus):
    
    A1 = np.linspace(0,1,len(fedata3[s:e+1])) 
    result50plus.append(np.interp(A, A1, fedata3[s:e+1]))
    resultbefore50.append(fedata3[s-5:s])
    resultafter50.append(fedata3[e:e+5])
resnp50plus = np.array(result50plus) 
mean50plus = np.nanmean(resnp50plus,axis = 0)



startbin75plus = [ x // 300 for x in plus75["ORF.start"]  ]
endbin75plus = [ x // 300 for x in plus75["ORF.end"] ]

A = np.linspace(0,1,100) 
result75plus = [] 
resultbefore75 = []
resultafter75 = []
for s,e in zip(startbin75plus,endbin75plus):
    
    A1 = np.linspace(0,1,len(fedata3[s:e+1])) 
    result75plus.append(np.interp(A, A1, fedata3[s:e+1]))
    resultbefore75.append(fedata3[s-5:s])
    resultafter75.append(fedata3[e:e+5])
resnp75plus = np.array(result75plus) 
mean75plus = np.nanmean(resnp75plus,axis = 0)



startbin100plus = [ x // 300 for x in plus100["ORF.start"]  ]
endbin100plus = [ x // 300 for x in plus100["ORF.end"] ]

A = np.linspace(0,1,100) 
result100plus = [] 
resultbefore100 = []
resultafter100 = []
for s,e in zip(startbin100plus,endbin100plus):
    
    A1 = np.linspace(0,1,len(fedata3[s:e+1])) 
    result100plus.append(np.interp(A, A1, fedata3[s:e+1]))
    resultbefore100.append(fedata3[s-5:s])
    resultafter100.append(fedata3[e:e+5])
resnp100plus = np.array(result100plus) 
mean100plus = np.nanmean(resnp100plus,axis = 0)



fig, axs = plt.subplots(2, 2, sharex = True, sharey = True)
axs[0, 0].plot(mean25plus)
fig.text(0.02, 0.5, 'Average FE usage, + strand', ha='center', va='center', rotation='vertical')
plt.xlabel("Gene bins - Chromosome 3", horizontalalignment='center', position = (0.01, 25))
plt.suptitle("Genes split into quantiles from lowest to highest RNA Levels")
axs[0, 0].set_title('0-25%')
axs[0, 1].plot(mean50plus, 'tab:orange')
axs[0, 1].set_title('25%-50%')
axs[1, 0].plot(mean75plus, 'tab:green')
axs[1, 0].set_title('50%-75%')
axs[1, 1].plot(mean100plus, 'tab:red')
axs[1, 1].set_title('75%-100%')
plt.show()


#minus strand quantiles

minus25 = q25chr3.loc[q25chr3['coding.strand']== '-']
minus50 = q50chr3.loc[q50chr3['coding.strand']== '-']
minus75 = q75chr3.loc[q75chr3['coding.strand']== '-']
minus100 = q100chr3.loc[q100chr3['coding.strand']== '-']




startbin25minus = [ x // 300 for x in minus25["ORF.start"]  ]
endbin25minus = [ x // 300 for x in minus25["ORF.end"] ]

A = np.linspace(0,1,100) 
result25 = [] 
resultbefore25 = []
resultafter25 = []
for s,e in zip(startbin25minus,endbin25minus):
    
    A1 = np.linspace(0,1,len(redata3[s:e+1])) 
    result25.append(np.interp(A, A1, redata3[s:e+1][::-1]))
    resultbefore25.append(redata3[s-5:s][::-1])
    resultafter25.append(redata3[e:e+5][::-1])
resnp25 = np.array(result25) 
mean25 = np.nanmean(resnp25,axis = 0)


startbin50minus = [ x // 300 for x in minus50["ORF.start"]  ]
endbin50minus = [ x // 300 for x in minus50["ORF.end"] ]

A = np.linspace(0,1,100) 
result50 = [] 
resultbefore50 = []
resultafter50 = []
for s,e in zip(startbin50minus,endbin50minus):
    
    A1 = np.linspace(0,1,len(redata3[s:e+1])) 
    result50.append(np.interp(A, A1, redata3[s:e+1][::-1]))
    resultbefore50.append(redata3[s-5:s][::-1])
    resultafter50.append(redata3[e:e+5][::-1])
resnp50 = np.array(result50) 
mean50 = np.nanmean(resnp50,axis = 0)



startbin75minus = [ x // 300 for x in minus75["ORF.start"]  ]
endbin75minus = [ x // 300 for x in minus75["ORF.end"] ]

A = np.linspace(0,1,100) 
result75 = [] 
resultbefore75 = []
resultafter75 = []
for s,e in zip(startbin75minus,endbin75minus):
    
    A1 = np.linspace(0,1,len(redata3[s:e+1])) 
    result75.append(np.interp(A, A1, redata3[s:e+1][::-1]))
    resultbefore75.append(redata3[s-5:s][::-1])
    resultafter75.append(redata3[e:e+5][::-1])
resnp75 = np.array(result75) 
mean75 = np.nanmean(resnp75,axis = 0)



startbin100minus = [ x // 300 for x in minus100["ORF.start"]  ]
endbin100minus = [ x // 300 for x in minus100["ORF.end"] ]

A = np.linspace(0,1,100) 
result100 = [] 
resultbefore100 = []
resultafter100 = []
for s,e in zip(startbin100minus,endbin100minus):
    
    A1 = np.linspace(0,1,len(redata3[s:e+1])) 
    result100.append(np.interp(A, A1, redata3[s:e+1][::-1]))
    resultbefore100.append(redata3[s-5:s][::-1])
    resultafter100.append(redata3[e:e+5][::-1])
resnp100 = np.array(result100) 
mean100 = np.nanmean(resnp100,axis = 0)


#reverse epsilon usage 

fig, axs = plt.subplots(2, 2, sharex = True, sharey = True)
axs[0, 0].plot(mean25)
fig.text(0.02, 0.5, 'Average RE usage, - strand', ha='center', va='center', rotation='vertical')
plt.xlabel("Gene bins - Chromosome 3", horizontalalignment='center', position = (0.01, 25))
plt.suptitle("Genes split into quantiles from lowest to highest RNA Levels")
axs[0, 0].set_title('0-25%')
axs[0, 1].plot(mean50, 'tab:orange')
axs[0, 1].set_title('25%-50%')
axs[1, 0].plot(mean75, 'tab:green')
axs[1, 0].set_title('50%-75%')
axs[1, 1].plot(mean100, 'tab:red')
axs[1, 1].set_title('75%-100%')
plt.show()

#######################################################################################################
startbins = [ x // 300 for x in chr1["ORF.start"]  ]
endbins = [ x // 300 for x in chr1["ORF.end"] ]



Z = np.linspace(0,1,100) 
result = [] 
for s,e in zip(startbins,endbins):
    Z1 = np.linspace(0,1,(e-s)+1) 
    result.append(np.interp(Z, Z1, fedata[s:e+1]))

resnp = np.array(result) 
meana = np.nanmean(resnp,axis = 0)



fig, axs = plt.subplots(1, squeeze=False)
axs[0,0].plot(meana)
plt.ylabel("Average Forward strand Epsilon usage")
plt.xlabel("Gene bins (mean)")
plt.show()




Z = np.linspace(0,1,100) 
result = [] 
for s,e in zip(startbins,endbins):
    Z1 = np.linspace(0,1,(e-s)+1) 
    result.append(np.interp(Z, Z1, redata[s:e+1]))

resnp = np.array(result) 
meana = np.nanmean(resnp,axis = 0)



fig, axs = plt.subplots(1, squeeze=False)
axs[0,0].plot(meana)
plt.ylabel("Average Reverse strand Epsilon usage")
plt.xlabel("Gene bins (mean)")
plt.show()








###############################################################################################################


startbins = [ x // 300 for x in chr1["ORF.start"]  ]
endbins = [ x // 300 for x in chr1["ORF.end"] ]



Z = np.linspace(0,1,100) 
result = [] 
for s,e in zip(startbins,endbins):
    Z1 = np.linspace(0,1,(e-s)+1) 
    result.append(np.interp(Z, Z1, cl_drt[s:e+1]))

resnp = np.array(result) 
meana = np.nanmean(resnp,axis = 0)



fig, axs = plt.subplots(1, squeeze=False)
axs[0,0].plot(meana)
plt.ylabel("Origin Efficiency")
plt.xlabel("Gene bins (mean)")
plt.show()




########################################################################################



Xbefore = np.linspace(-0.1,0,6)
Xafter = np.linspace(1,1.1,6)

Result = []
resultbefore = []
resultafter = []
for s,e in zip(startbinplus,endbinplus):
    if e-s == 0:
     continue
    Z1 = np.linspace(0,1,(e-s)+1)  
    # fedata for + strand
    Result.append(np.interp(Z, Z1, fedata[s:e+1]))
    resultbefore.append( fedata[s-5:s+1] )  #5 bins before start position
    resultafter.append( fedata[e:e+5+1]) #and the bins after. 
   
resnp1 = np.array(Result)
mean1f = np.nanmean(resnp1, axis = 0)

resnp2 = np.array(resultbefore)
mean2 = np.nanmean(resnp2, axis = 0)

resnp3 = np.array(resultafter)
mean3 = np.nanmean(resnp3, axis = 0)
    


fig, (ax20, ax21) = plt.subplots(nrows=2, ncols=1, sharex = True)
ax20.plot(Xbefore, mean2)
ax21.plot(Xafter, mean3)
plt.ylabel("Forward Epsilon usage(plus strand)")
plt.xlabel("Gene bins")
plt.title("1.5kb before/after the position of genes ")
plt.show()

fig,(ax22) = plt.subplots(nrows = 1, ncols = 1)
ax22.plot(Z,mean1f)
plt.ylabel("Average forward epsilon usage")
plt.xlabel("Gene Bins (mean)")



Result1 = []
resultbefore1 = []
resultafter1 = []
for s,e in zip(startbinminus,endbinminus):
    if e-s == 0:
     continue
    Z1 = np.linspace(0,1,(e-s)+1) 
    Result1.append(np.interp(Z, Z1, redata[s:e+1][::-1])) #in reverse order, so the ::-1.
    resultbefore1.append( redata[e:e+5+1][::-1]  )  #Here I add the 5 bins before start position
    resultafter1.append( redata[s-5:s+1][::-1])  #and the bins after.





resnp4 = np.array(resultbefore1[:-1]) 

mean4 = np.nanmean(resnp4,axis = 0)


resnp5 = np.array(resultafter1[1:]) 
mean5 = np.nanmean(resnp5,axis = 0)


resnp6 = np.array(Result1) 
mean6 = np.nanmean(resnp6,axis = 0)


fig, (ax23, ax24) = plt.subplots(nrows=2, ncols=1, sharex = True)
ax23.plot(Xbefore, mean4)
plt.ylabel("Reverse Epsilon Usage")
ax24.plot(Xafter, mean5)
plt.ylabel("Reverse Epsilon usage(minus strand)")
plt.xlabel("Genes bins")
plt.show()



plt.plot(Xbefore,mean4, 'b')
plt.plot(Xafter,mean5, 'b')
plt.plot(Z, mean6, 'b')
plt.ylabel("Reverse Epsilon usage(minus strand)")
plt.xlabel("Genes bins")
plt.title("1.5kb before/after the position of genes ")
plt.show()

plt.plot(Xbefore,mean2, 'b')
plt.plot(Xafter,mean3, 'b')
plt.plot(Z, mean1f, 'b')
plt.ylabel("Forward Epsilon usage(plus strand)")
plt.xlabel("Gene bins")
plt.title("1.5kb before/after the position of genes ")
plt.show()


fig,(ax25) = plt.subplots(nrows=1,ncols=1)
plt.plot(Z, mean1f)
#origin efficiencies - upstream/downstream of genes?

fig,(ax26) = plt.subplots(nrows=1,ncols=1)
plt.plot(Z, mean6)
plt.ylabel("Reverse Epsilon usage(plus strand)")
plt.xlabel("Gene bins")
##########################################################################################################



Result = []
resultbefore = []
resultafter = []
for s,e in zip(startbinplus,endbinplus):
    if e-s == 0:
     continue
    Z1 = np.linspace(0,1,(e-s)+1)  
    # fedata for + strand
    Result.append(np.interp(Z, Z1, cl_drt[s:e+1]))
    resultbefore.append( cl_drt[s-6:s+1] )  #5 bins before start position
    resultafter.append( cl_drt[e:e+6+1]) #and the bins after. 
   
resnp1 = np.array(Result)
mean1 = np.nanmean(resnp1, axis = 0)

resnp2 = np.array(resultbefore)
mean2 = np.nanmean(resnp2, axis = 0)

resnp3 = np.array(resultafter)
mean3 = np.nanmean(resnp3, axis = 0)
    


fig, (ax25, ax26) = plt.subplots(nrows=2, ncols=1, sharex = True)
ax25.plot(Xbefore, mean2)
ax26.plot(Xafter, mean3)
plt.ylabel("Efficiency")
plt.xlabel("Gene bins")
plt.show()


plt.plot(Xbefore, mean2 , 'b')
plt.plot(Xafter, mean3, 'b')
plt.plot(Z, meana,'b')
plt.xlabel("Gene Bins")
plt.ylabel("Origin Efficiency")
plt.show()



















 