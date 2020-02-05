# original code from Wu Tongfei
# https://github.com/Ngtungfei/Graphene-Oxide-Generator2.git
#
# small amendments by V Erastova, University of Edinburgh, 5th Feb 2020
#

import matplotlib.pyplot as plt
import numpy as np
import os
from collections import Counter
import random
from matplotlib.lines import Line2D
from tkinter import filedialog
#import tkinter.messagebox
from tkinter import *
from Gopackages.hexahapto  import buildhexahapto as bup
from Gopackages.graphene import buildgraphene as bug
#from gopackages.parallelogram import buiildarallelogram as bupar


class Example(Frame):

    def __init__(self, parent):
        Frame.__init__(self, parent)   
        self.parent = parent
        self.initUI()
        
    def initUI(self):
        self.parent.title("GO")
        sw=self.parent.winfo_screenwidth()
        sh=self.parent.winfo_screenheight()
        x = int(self.parent.winfo_screenwidth()*0.8/5*4)
        y = (self.parent.winfo_screenheight()-720)/2
        self.parent.geometry('+%d+%d' % (x, y))
        self.pack(fill=BOTH, expand=1)
        
        #----------------bong length /Angstroms
        self.cabond=1.34
        self.ccbond= 1.54
        self.cdobond=1.23
        self.cobond=1.43
        self.ohbond=0.96

        self.wb = Label(self, text='Graphene Size: ',width=18,anchor=E)
        self.wb.grid(row=0, column=0)        
        self.wb = Label(self, text='*compulsory actions ',width=18,anchor=E, fg='red')
        self.wb.grid(row=0, column=1)      
        self.wb = Label(self, text='x (nm)',width=18,anchor=E)
        self.wb.grid(row=1, column=0)
        self.wb = Label(self, text='y (nm)',width=18,anchor=E)
        self.wb.grid(row=2, column=0)

        self.ww = Scale(self, from_=1, to=11, resolution=0.5, orient=HORIZONTAL, fg='red')
        self.ww.grid(row=1, column=1)
        self.www = Scale(self,from_=1, to=11,resolution=0.5, orient=HORIZONTAL, fg='red')
        self.www.grid(row=2, column=1)
        
        self.check1=IntVar()
        self.che1 = Checkbutton(self, text="Hexagonal (Radii=x) ",variable=self.check1)
        self.che1.grid(row=3, column=0)
        
        self.quitButton1 = Button(self, text = '*1.Generate',width=15, fg='red', command = self.newwindows)
        self.quitButton1.grid(row=3, column=1)
      

        self.wb = Label(self, text='Edge Carbons (%)',width=18,anchor=E)
        self.wb.grid(row=5, column=0)
        self.wb = Label(self, text='Bulk Carbons (%)',width=18,anchor=E)
        self.wb.grid(row=6, column=0)
        self.wb = Label(self, text='Total rate (%)',width=18,anchor=E)
        self.wb.grid(row=7, column=0)

        self.w1 = Scale(self, from_=0, to=100,orient=HORIZONTAL, fg='red')
        self.w1.grid(row=5, column=1)
        self.w2 = Scale(self, from_=0, to=100,orient=HORIZONTAL, fg='red')
        self.w2.grid(row=6, column=1)
        self.w3 = Scale(self, from_=0, to=100,orient=HORIZONTAL, fg='red')
        self.w3.grid(row=7, column=1)

        self.quitButton2 = Button(self, text = '*2.Delete Carbons',width=15, fg='red',command = self.remcarbons)
        self.quitButton2.grid(row=8, column=1)
       
        self.wc = Label(self, text='Add to Edge C:',width=18,anchor=E)
        self.wc.grid(row=9, column=0)
        
        self.quitButton3 = Button(self, text = '3.-COOH and -C=O',width=15,fg='blue', command = self.addco)
        self.quitButton3.grid(row=9, column=1)

        self.wd = Label(self, text='Add to C=C bonds:',width=18,anchor=E)
        self.wd.grid(row=10, column=0)
        
        self.w6 = Scale(self, from_=0, to=100,orient=HORIZONTAL, fg='blue')
        self.w6.grid(row=10, column=1)
        
        self.quitButton3 = Button(self, text = '4. -O- (%)',width=15,fg='blue', command = self.addo)
        self.quitButton3.grid(row=11, column=1)

        self.w5 = Scale(self, from_=0, to=100,orient=HORIZONTAL, fg='blue')
        self.w5.grid(row=12, column=1)
        
        self.quitButton3 = Button(self, text = '5. -OH (%)',width=15,fg='blue', command = self.addoh)
        self.quitButton3.grid(row=13, column=1)

        self.check=IntVar()
        self.che = Checkbutton(self, text="mol2 force field",variable=self.check)
        self.che.grid(row=14, column=0)

        self.b1=Button(self, text="SAVE .mol2",width=15,command= self.save)
        self.b1.grid(row=15, column=1)

    def newwindows(self):
        self.xmin=0
        self.ymin=0
        self.xmax=0
        self.ymax=0

        self.cno=0
        self.clo=0
        self.cono=0
        self.coohno=0
        self.ohno=0
        self.ono=0
        self.step=1
        
        if self.check1.get():
            n=((int(self.ww.get()*10)//(1.732*self.cabond))//2)*2+1
            self.npatomsdata,self.npbondsdata,self.xmin,self.xmax,self.ymin,self.ymax=bup(n)
            self.cno=self.npatomsdata.shape[0]
            self.clo= self.cno
            self.plot(self.npatomsdata,self.npbondsdata)
        
        else:
            n=((int((self.www.get()*10)//(0.866*self.cabond))-3)//2)*2+3
            m=int((self.ww.get()*10)//(3*self.cabond)//2)*2*2
            self.npatomsdata,self.npbondsdata,self.xmin,self.xmax,self.ymin,self.ymax=bug(n,m)
            self.cno=self.npatomsdata.shape[0]
            self.clo= self.cno
            self.plot(self.npatomsdata,self.npbondsdata)

    def plot(self,npatomsdata,npbondsdata):
        try:
            self.fig.clf()            
            self.fig = plt.figure()
        except:
            self.fig = plt.figure()
            #self.fig, ax = plt.subplots(1, 1)
            
        ax = self.fig.add_subplot(111)
        ax.set_xlim([self.xmin-5,self.xmax+5])
        ax.set_ylim([self.ymin-5,self.ymax+5])
        
        dic={1:'c',2:'k',3:'r',4:'b'}
        col={'C':'k','O':'r','H':'b'}
        siz={'C':'600','O':'500','H':'300'}
        
        patchC = Line2D(range(1), range(1), color="white", marker='o', markersize=10,markeredgecolor="black", markerfacecolor="black",alpha=0.6)
        patchO = Line2D(range(1), range(1), color="white", marker='o', markersize=9,markeredgecolor="red", markerfacecolor="red",alpha=0.6)
        patchH = Line2D(range(1), range(1), color="white", marker='o', markersize=6, markeredgecolor="blue",markerfacecolor="blue",alpha=0.3) 
        plt.legend([patchC,patchO,patchH],['C','O','H'],loc='lower left',bbox_to_anchor=(1, 0),numpoints = 1)
        a=self.ww.get()
        b=self.www.get()
        
        if self.check1.get():
            b=self.ww.get()
        for i in npbondsdata:
            ax.plot([float(npatomsdata[i[0]-1][2]),float(npatomsdata[i[1]-1][2])],[float(npatomsdata[i[0]-1][3]),float(npatomsdata[i[1]-1][3])],color=dic[i[2]])
            ax.scatter(float(npatomsdata[i[0]-1][2]),float(npatomsdata[i[0]-1][3]),color=col[npatomsdata[i[0]-1][1]],s=int(siz[npatomsdata[i[0]-1][1]])/a/b/1.5,alpha=.4)
            ax.scatter(float(npatomsdata[i[1]-1][2]),float(npatomsdata[i[1]-1][3]),color=col[npatomsdata[i[1]-1][1]],alpha=.3,s=int(siz[npatomsdata[i[1]-1][1]])/a/b/1.5)
        
        ax.set_aspect('equal')
        ax.set_xlabel('Angstroms')
        ax.set_ylabel('Angstroms')
        ax.set_title('%d Carbons %d left:\n %d(=O) %d(COOH) %d(-O-) %d(OH)'%(self.cno,self.clo,self.cono,self.coohno,self.ono,self.ohno))
        plt.show()

    # decide on edges and bulk C
    def edgebulk(self,npatomsdata,npbondsdata):
        count0=Counter(npbondsdata[:,0])
        count1=Counter(npbondsdata[:,1])
        removeatoms=[]
        edgeatoms=[]
        bulkatoms=[]
        for i in range(npatomsdata.shape[0]):
            count=count0[i+1]+count1[i+1]
            if count ==1:
                removeatoms.append(i+1)
            elif count ==2:
                edgeatoms.append(i+1)
            elif count ==3:
                bulkatoms.append(i+1)
        return removeatoms,edgeatoms,bulkatoms

    def edges(self,edgeatoms):
        edgeatomco=[]
        edgeatomcooh=[]
        for i in edgeatoms:
            k=0
            for j in self.npbondsdata2:
                if i ==j[0] or i ==j[1]:
                    k+=j[2]
            if k==2:
                edgeatomco.append(i)
            elif k==3:
                edgeatomcooh.append(i)
        return edgeatomcooh,edgeatomco

    # remove C to make holes
    def remcarbons(self):
        self.step=2
        s=0
        bonds=self.npbondsdata.tolist()
        self.npbondsdata2=np.array(bonds)

        
        if self.w1.get()+self.w2.get()==0:
        	#if no holes and edges to be made 
            print('no edges, no bulk atoms affected')        
            print('nothing to do')
            return

                        
        if self.w3.get()==0:
        	#if no holes and edges to be made 
            print('will only make a small single hole')
            print('better set the rate above zero')
            return

        while 1:
            removeatoms,edgeatoms,bulkatoms=self.edgebulk(self.npatomsdata,self.npbondsdata2)
            bins=[]           
            remc=np.random.choice(edgeatoms, 1)
            if random.randrange(0,100) < self.w1.get():
                for j in self.npbondsdata2:
                    if remc == j[0] or remc == j[1]:
                        bins.append(j.tolist())
            bonds=self.npbondsdata2.tolist()
            for i in bins:
                for j in bonds:
                    if i==j:
                        bonds.remove(j)
            self.npbondsdata2=np.array(bonds)
            
            removeatoms,edgeatoms,bulkatoms=self.edgebulk(self.npatomsdata,self.npbondsdata2)
            self.removelone(removeatoms)
            s=1-self.npbondsdata2.shape[0]/self.npbondsdata.shape[0]
            
            if s>self.w3.get()/100:
                break

            removeatoms,edgeatoms,bulkatoms=self.edgebulk(self.npatomsdata,self.npbondsdata2)
            bins=[]           
            remc=np.random.choice(bulkatoms, 1)
            if random.randrange(0,100) < self.w2.get():
                for j in self.npbondsdata2:
                    if  remc == j[0] or remc == j[1]:
                        bins.append(j.tolist())
            bonds=self.npbondsdata2.tolist()
            for i in bins:
                for j in bonds:
                    if i==j:
                        bonds.remove(j)
            self.npbondsdata2=np.array(bonds)
            removeatoms, edgeatoms,bulkatoms=self.edgebulk(self.npatomsdata,self.npbondsdata2)
            self.removelone(removeatoms)
            s=1-self.npbondsdata2.shape[0]/self.npbondsdata.shape[0]
            
            if s>self.w3.get()/100:
                break
           
        removeatoms, edgeatoms,bulkatoms=self.edgebulk(self.npatomsdata,self.npbondsdata2)
        while removeatoms!=[]:
            self.removelone(removeatoms)
            removeatoms,edgeatoms,bulkatoms=self.edgebulk(self.npatomsdata,self.npbondsdata2)
            
        self.clo=len(edgeatoms)+len(bulkatoms)
        print('removed Cs and some bonds between')

        self.plot(self.npatomsdata,self.npbondsdata2)

    def removelone(self,removeatoms):
        bins=[]
        for i in removeatoms:
            for j in self.npbondsdata2:
                if i == j[0] or i == j[1]:
                    bins.append(j.tolist())    
        bonds=self.npbondsdata2.tolist()
        for i in bins:
            for j in bonds:
                if i==j:
                    bonds.remove(j)
        self.npbondsdata2=np.array(bonds)
                    

    # add terminations onto holes and edges    
    def addco(self):
        self.step=3
        removeatoms,edgeatoms,bulkatoms=self.edgebulk(self.npatomsdata,self.npbondsdata2)
        edgeatomcooh,edgeatomco=self.edges(edgeatoms)
        
        npatomsdatas=self.npatomsdata.tolist()
        npbondsdata2s=self.npbondsdata2.tolist()
        self.cono=len(edgeatomco)

        for i in edgeatomco:
            npatomsdatas[i][5]='CG2O5'
            npatomsdatas[i][8]='0.47'
            a=[]
            for j in self.npbondsdata2:
                if i ==j[0] :
                    a.append(int(j[1]))
                elif i ==j[1]:
                    a.append(int(j[0]))
                    
            bx=(float(self.npatomsdata[a[0]-1][2])+float(self.npatomsdata[a[1]-1][2]))/2
            by=(float(self.npatomsdata[a[0]-1][3])+float(self.npatomsdata[a[1]-1][3]))/2
            cx=float(self.npatomsdata[i-1][2])
            cy=float(self.npatomsdata[i-1][3])
            
            x=cx-2*(bx-cx)*self.cdobond/self.cabond
            y=cy-2*(by-cy)*self.cdobond/self.cabond
            
            taotal=len(npatomsdatas)
            npbondsdata2s.append([i, taotal+1, 3])
            npatomsdatas.append([taotal+1,'O', x, y, '0.0', 'OG2D3', '1', 'GO', '-0.47'])
        self.npatomsdata=np.array(npatomsdatas)
        self.npbondsdata3=np.array(npbondsdata2s)
        
        print('added C=O')

        self.addcoohs(edgeatomcooh)
      
    #add terminations to edges and holes
    def addcoohs(self,edgeatomcooh):
        npatomsdatas=self.npatomsdata.tolist()
        npbondsdata2s=self.npbondsdata3.tolist()
        self.coohno=len(edgeatomcooh)
        for i in edgeatomcooh:
            a=[]
            for j in self.npbondsdata2:
                if i ==j[0] :
                    a.append(int(j[1]))
                elif i ==j[1]:
                    a.append(int(j[0]))
                    
            bx=(float(self.npatomsdata[a[0]-1][2])+float(self.npatomsdata[a[1]-1][2]))/2
            by=(float(self.npatomsdata[a[0]-1][3])+float(self.npatomsdata[a[1]-1][3]))/2
            cx=float(self.npatomsdata[i-1][2])
            cy=float(self.npatomsdata[i-1][3])
            x=cx-2*(bx-cx)*self.ccbond/self.cabond
            y=cy-2*(by-cy)*self.ccbond/self.cabond
            c=random.sample([-1.0,1.0], 1)
            c=c[0]
            
            taotal=len(npatomsdatas)
            npbondsdata2s.append([i, taotal+1, 1])
            npatomsdatas.append([taotal+1,'C', x, y, '0.0', 'CG2O2', '1', 'GO', '0.46'])

            x=cx-2*(bx-cx)*(self.ccbond+self.cdobond/2)/self.cabond
            y=cy-2*(by-cy)*(self.ccbond+self.cdobond/2)/self.cabond

            taotal=len(npatomsdatas)
            npbondsdata2s.append([taotal, taotal+1, 3])
            npatomsdatas.append([taotal+1,'O', x, y, -1*c*self.cdobond*0.866, 'OG2D1', '1', 'GO', '-0.46'])

            x=cx-2*(bx-cx)*(self.ccbond+self.cobond/2)/self.cabond
            y=cy-2*(by-cy)*(self.ccbond+self.cobond/2)/self.cabond

            taotal=len(npatomsdatas)
            npbondsdata2s.append([taotal-1, taotal+1, 1])
            npatomsdatas.append([taotal+1,'O', x, y, c*self.cobond*0.866, 'OG311', '1', 'GO', '-0.51'])

            x=cx-2*(bx-cx)*(self.ccbond+self.cobond/2+self.ohbond)/self.cabond
            y=cy-2*(by-cy)*(self.ccbond+self.cobond/2+self.ohbond)/self.cabond

            taotal=len(npatomsdatas)
            npbondsdata2s.append([taotal, taotal+1, 4])
            npatomsdatas.append([taotal+1,'H', x, y, c*self.cobond*0.866, 'HGP1', '1', 'GO', '0.51'])
            
        self.npatomsdata=np.array(npatomsdatas)
        self.npbondsdata3=np.array(npbondsdata2s)
        
        print('and added C(=O)OH')
        self.plot(self.npatomsdata,self.npbondsdata3)
        


    # add bridging -O- to C=C 
    def addo(self):
        if self.w6.get()==0:
            return

        #if C=O COOH were not added to edges
        if not hasattr( self,'npbondsdata3'):
        	print('you have skipped step 3 - not added carboxilates to holes ')
        	self.npbondsdata3=self.npbondsdata2


        self.step=4
        npatomsdatas=self.npatomsdata.tolist()
        npbondsdata2s=self.npbondsdata3.tolist()
        doublebonds=[]
        for i in range(len(npbondsdata2s)):
            if npbondsdata2s[i][2]==2:
                doublebonds.append(i)
        n=int(self.w6.get()*len(doublebonds)/100)
        remc=random.sample(doublebonds, n)
        self.ono=len(remc)
        for i in remc:
            a=npbondsdata2s[i][0]
            b=npbondsdata2s[i][1]

            npatomsdatas[a-1][5]='CG3C31'
            npatomsdatas[a-1][8]='0.20'
            npatomsdatas[b-1][5]='CG3C31'
            npatomsdatas[b-1][8]='0.20'

            npbondsdata2s[i][2]=1

            x=(float(self.npatomsdata[a-1][2])+float(self.npatomsdata[b-1][2]))/2
            y=(float(self.npatomsdata[a-1][3])+float(self.npatomsdata[b-1][3]))/2
            
            taotal=len(npatomsdatas)
            npbondsdata2s.append([a, taotal+1, 3])
            npbondsdata2s.append([b, taotal+1, 3])
            c=random.sample([-1.0,1.0], 1)
            c=c[0]
            npatomsdatas.append([taotal+1,'O', x, y, str(c), 'OG3C31', '1', 'GO', '-0.40'])
            
        self.npatomsdata=np.array(npatomsdatas)
        self.npbondsdata4=np.array(npbondsdata2s)
        print('added' ,n ,' -O-')

        self.plot(self.npatomsdata,self.npbondsdata4)
    
    # add -OH to C=C
    def addoh(self):
        if self.w5.get()==0:
            return
        if not hasattr( self,'npbondsdata3'):
        	print('you have skipped step 3 - not added carboxilates to holes ')
        	self.npbondsdata3=self.npbondsdata2

        if not hasattr( self,'npbondsdata4'):
        	print('you have skipped step 4 - not added bridging C ')
        	self.npbondsdata4=self.npbondsdata3


        self.step=5
        npatomsdatas=self.npatomsdata.tolist()
        npbondsdata2s=self.npbondsdata4.tolist()
        doublebonds=[]
        for i in range(len(npbondsdata2s)):
            if npbondsdata2s[i][2]==2:
                doublebonds.append(i)
        n=int(self.w5.get()*len(doublebonds)/100)
        remc=random.sample(doublebonds, n)
        self.ohno=len(remc)*2
        for i in remc:
            a=npbondsdata2s[i][0]
            b=npbondsdata2s[i][1]

            npatomsdatas[a-1][5]='CG3C51'
            npatomsdatas[a-1][8]='0.14'
            npatomsdatas[b-1][5]='CG3C51'
            npatomsdatas[b-1][8]='0.14'

            npbondsdata2s[i][2]=1

            x=(float(self.npatomsdata[a-1][2])+float(self.npatomsdata[b-1][2]))/2
            y=(float(self.npatomsdata[a-1][3])+float(self.npatomsdata[b-1][3]))/2
            
            c=random.sample([-1.0,1.0], 1)
            c=c[0]
            taotal=len(npatomsdatas)
            npbondsdata2s.append([a, taotal+1, 4])
            npatomsdatas.append([taotal+1,'O', x, y, str(c), 'OG311', '1', 'GO', '-0.65'])

            taotal=len(npatomsdatas)
            npbondsdata2s.append([taotal,  taotal+1, 4])
            npatomsdatas.append([taotal+1,'H', x, y, str(c*2), 'HGP1', '1', 'GO', '0.65'])

            taotal=len(npatomsdatas)
            npbondsdata2s.append([b, taotal+1, 4])
            npatomsdatas.append([taotal+1,'O', x, y, str(-1*c), 'OG311', '1', 'GO', '-0.65'])

            taotal=len(npatomsdatas)
            npbondsdata2s.append([taotal, taotal+1, 4])
            npatomsdatas.append([taotal+1,'H', x, y, str(-1*c*2), 'HGP1', '1', 'GO', '0.65'])

        self.npatomsdata=np.array(npatomsdatas)
        self.npbondsdata5=np.array(npbondsdata2s)
        print('added ', n ,' -OH')

        self.plot(self.npatomsdata,self.npbondsdata5)

    def save(self):
        typedic={'CG2O5':'C.2','CG2O2':'C.2','OG2D1':'O.2','OG311':'O.3','HGP1':'H','OG2D3':'O.2','CG2R61':'C.2','CG3C51':'C.3','CG3C31':'C.3','OG3C31':'O.3'}
        npatomsdatas=self.npatomsdata.tolist()
        if self.step ==1:
            npbondsdata2s=self.npbondsdata.tolist()        
        elif self.step ==2:
            npbondsdata2s=self.npbondsdata2.tolist()
        elif self.step ==3:
            npbondsdata2s=self.npbondsdata3.tolist()
        elif self.step ==4:
            npbondsdata2s=self.npbondsdata4.tolist()
        elif self.step ==5:
            npbondsdata2s=self.npbondsdata5.tolist()
                  
        ftypes=[('Tripo mol2 file','.mol2')]
        f = filedialog.asksaveasfilename(filetypes=ftypes,initialdir=os.curdir,initialfile='GO.mol2')
        if f=='':
            return
        file = open(f, "w")
        k=[]
        for i in range(len(npatomsdatas)):
            for j in npbondsdata2s:
                if i+1==j[0] or i+1 ==j[1]:
                    k.append(i+1)
        count=Counter(np.array(k))
                    
        file.write('@<TRIPOS>MOLECULE\n')
        file.write('Graphene Oxide # %d Carbon %d left: %d(=O) %d(COOH) %d(-O-) %d(-OH)\n'%(self.cno,self.clo,self.cono,self.coohno,self.ono,self.ohno))
        file.write('%d %d 0 0 0\n'%(len(count.keys()),len(npbondsdata2s)))
        file.write('SMALL\n')
        file.write('GASTEIGER\n')
        file.write('\n')
        file.write('@<TRIPOS>ATOM\n')
        opp=0
        aDict = {}

        for i in count.keys():
            j = npatomsdatas[i-1]
            opp+=1
            if self.check.get():
                file.write('%d %s %.4f %.4f %.4f %s 1 GO %.4f\n'%(opp,j[1],float(j[2]),float(j[3]),float(j[4]),typedic[j[5]],float(j[8])))
            else:
                file.write('%d %s %.4f %.4f %.4f %s 1 GO %.4f\n'%(opp,j[1],float(j[2]),float(j[3]),float(j[4]),j[5],float(j[8])))
            aDict[str(j[0])] = opp
                    
        file.write('\n@<TRIPOS>BOND\n')
        for i in range(len(npbondsdata2s)):
            if npbondsdata2s[i][2]==2 or npbondsdata2s[i][2]==3:
                file.write('%d %d %d %d\n'%(i+1,aDict[str(npbondsdata2s[i][0])],aDict[str(npbondsdata2s[i][1])],2))
            else:
                file.write('%d %d %d %d\n'%(i+1,aDict[str(npbondsdata2s[i][0])],aDict[str(npbondsdata2s[i][1])],1))
        file.close()

def main():
  
    root = Tk()
    ex = Example(root)
    root.mainloop()
    
if __name__ == '__main__':
    main()  



    
    


        
 

