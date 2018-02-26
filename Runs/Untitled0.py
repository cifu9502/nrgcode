
# coding: utf-8

# In[1]:

from __future__ import division 
#get_ipython().magic(u'pylab inline')
import matplotlib
import numpy as np
import string
import os
txt = 'run'
txt2 = '/rho_'
densPath = '/home/cifucito/nrgcode/TwoChNRG/Images/Density/run'


# In[2]:

for j in range(6,7):
    for i in range(4):
        output = txt+repr(j)+txt2 +repr(i)+'_'+repr(i)+'_OmegaRhow.dat'
        #output = '/home/cifucito/nrgcode/TwoChNRG/src/Main/Run/Run2DOtM/rho_2_2_OmegaRhow.dat'
        print output
        infile = open(os.path.abspath(output), 'r')
        text = infile.readlines()
        vec = []
        vec2 = []
        for x in text:
            #print(list(x))
            a =x.split(' ')
            vec.append(float(a[-3]))
            vec2.append(float(a[0]))
        vec = np.array(vec)    
        maxima = vec[r_[True, vec[1:] > vec[:-1]] & r_[vec[:-1] > vec[1:], True]]
        maxima = np.sort(maxima)
        plt.plot(vec2,vec)    
        plt.xlabel('$\omega$',fontsize=19)
        plt.ylabel('Spectral Density ',fontsize=18)
        plt.title(" , ".join(str(x) for x in maxima[-3:]),fontsize=10)
        directory = densPath + repr(j)
        #get_ipython().system(u'mkdir {directory}')
        save = directory +'/'+txt2 +repr(i)+'_'+repr(i)+'_OmegaRhow.png'
       
        
        plt.savefig(save)
        close()


# In[7]:

#get_ipython().system(u'mkdir {txt}')


# In[18]:

#vec


# In[41]:




# In[41]:




# In[ ]:



