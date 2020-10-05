#################################################################################################################
###################################引入模块区域##################################################################
import xlrd#读Excle
import numpy as np
import re###正则表达式
###分割线
import random
import matplotlib.pyplot as plt #绘制图像要用到的pyplot
import itertools
import math
import copy
##取消科学计数法https://blog.csdn.net/Andrew_jdw/article/details/82350041
np.set_printoptions(suppress=True)
###################################引入模块区域##################################################################
#################################################################################################################


#####第一步
####读写Excel文件
#################################################################################################################
#################全局变量：开始#######################################
#################全局变量：开始#######################################
######相当于一个登记区域

####引用：https://www.runoob.com/python/python-reg-expressions.html
####引用：https://blog.csdn.net/u010412858/article/details/83062200
wb = xlrd.open_workbook("更新表格信息.xlsx")#打开文件
##print(wb.sheet_names())#获取所有表格名字
sheet1 = wb.sheet_by_index(1)#通过索引获取表格

############################################
###处理不同表格都需要用的变量###
inputstr=[]
nrow=sheet1.nrows
##print(nrow)
haveRecord=False
####汇总数据变量
Zxqdnum=0###总需求点数量
Zxqdzl=0###总需求点重量
Zxqdtj=0###总需求点体积
Zclnum=[]###总车辆个数
Zclzl=[]###总车辆重量
Zcltj=[]###总车辆体积
Zgylzl=0###供应量重量
Zgyltj=0###供应量体积 
############################################

############################################
###处理需求点1到需求点n的分表需要用的变量###
###装需求点的矩阵,按顺序
SPTSACDarr=[]#########################################重要
pattern = re.compile(r' *需求点\d')
numpattn=re.compile('\d')
sptsarr=[]###装需求点的矩阵,顺序是打乱的
spots=[]###需求点顺序表，与矩阵一一对应
############################################

############################################
###处理车辆类型要用到的变量###
patternCLLX = re.compile(r' *车辆类型')
###车辆类型表
CLLXarr=[]############################################重要
CLLXtypenum=[]########################################重要
CLLXtitl=[]###########################################重要
##############################################

############################################
###处理距离矩阵要用到的变量###
patternJL = re.compile(r'距离')
###距离表
JLarr=[]##############################################重要
############################################
###处理供应量矩阵要用到的变量###
patternGYL=re.compile(r'供应量')
###供应量表
GYLarr=[]#############################################重要

############################################
#####总表开始
#####需求点汇总
XQDHZarr=[]###########################################重要
#####物资类型汇总
WZLXHZarr=[]##########################################重要
#####总表结束
############################################

#################全局变量：结束#######################################
#################全局变量：结束#######################################


##########========除掉List中某些元素：该元素只含有""的除空函数开始
def chukong(Alist):
    Blist= [i for i in Alist if set(i)!={""}]
    return np.array(Blist)
##########========除掉List中某些元素：该元素只含有""的除空函数结束

#########============处理记录函数开始
def chuliJilu():
    ######申明全局变量，可以在函数外面改变变量的值########
    global pattern,patternCLLX,patternJL
    global haveRecord,inputstr,spots,sptsarr,CLLXarr,JLarr
    global CLLXtypenum,CLLXtitl
    global patternGYL,GYLarr
        
    ###处理类型标志初始化#####
    signXQD=False ###处理需求点1到需求点n的标志
    signCLLX=False ###处理车辆类型的标志
    signJL=False ###处理距离矩阵的标志   
    signGYL=False###处理供应量的标志
    ###切分数据标志初始化####
    sep=0###分隔列的位置
    start=[]###开始位置
    findblank=False###未找到第一个有效空格
##    print("处理记录")
    ##==
    arrrow=len(inputstr)####有多少行
##    print("行数：",arrrow)
    inputarray=np.array(inputstr)
    firstrow=inputarray[0,:]
    
##    print("检查第一行，也就是标题行")
##    print("firstrow",firstrow)
    try:
        for idx,i in enumerate(firstrow):###获取index的同时获取数值,参考：https://www.cnblogs.com/nlyangtong/p/11905028.html
                ###以下代码在标题行（第一行）中匹配关键字确定表格的类型
                Obj=pattern.findall(i)###单个需求点                     
                objCLLX=patternCLLX.findall(i)###车辆类型
                objJL=patternJL.findall(i)###距离                
                objGYL=patternGYL.findall(i)###供应量
                if Obj !=[]:
##                    print(Obj)
                    signXQD=True
                    spots.append(Obj[0])###将需求点加入到顺序表中
                    start.append(idx)###将开始位置记录下来
                if objCLLX !=[]:
##                    print(objCLLX)
                    signCLLX=True
                    start.append(idx)###将开始位置记录下来
                if objJL !=[]:
##                    print(objJL)
                    signJL=True
                    start.append(idx)###将开始位置记录下来
                if objGYL !=[]:
##                    print(objGYL)
                    signGYL=True
                    start.append(idx)
                
    except:##遇到空，跳过
        pass
    ####下面读取start的数值
##    print(start)

    ####下面判断分隔符sep的位置，以""作为分隔符
    for idx, i in enumerate(firstrow):
        if idx > start[0] and i=="" and  findblank==False :###不是开始的空格符,而且还未找到有效空格
            sep=idx
            findblank=True
##    print("分隔符位置：",sep)
    ####下面获取数据长度
    datalen=sep-1-start[0]###数据长度
##    print("数据长度",datalen)
    
    ###下面是根据对表格类型的判断结果，分类别进行处理
##    print("分类别处理")
    if signXQD==True:
##        print("处理单个需求点数据")
        stlen=len(start)##一行有几个需求点的数据
        sti=0
        while sti < stlen:
            spot=inputarray[1:arrrow,start[sti]+1:start[sti]+1+datalen]###需求点数据
            spot=np.array(spot)
            hi=spot.shape[0]##行
            li=spot.shape[1]##列
            spotdata=np.zeros((hi,li),dtype=float)
##            print(hi,li)
            ihi=0            
            while ihi<hi:
                ili=0
                while ili< li:                    
                    if spot[ihi,ili]=="":###把空的换成"0"
                        spotdata[ihi,ili]=int(0)##矩阵中只要有小数，整数也会转换成小数
                    else:
                        spotdata[ihi,ili]=float(spot[ihi,ili])
                    ili +=1
                ihi +=1
            sptsarr.append(spotdata)###把需求点数据矩阵保存好
            sti +=1
##        print("处理完毕")
    elif signCLLX==True:
        CLLXtypenum.append(cllxtypenum(inputarray,start,datalen))
        CLLXtitl.append(cllxtitle(inputarray,start,datalen))
##        print(CLLXtypenum)
##        print("处理车辆类型数据")
        spC=np.array(inputarray[1:arrrow,start[0]+1:start[0]+1+datalen])
        hi=spC.shape[0]##行
        li=spC.shape[1]##列
        spCdata=np.zeros((hi,li))
##        print(hi,li)
        ihi=0        
        while ihi<hi:
            ili=0
            while ili< li:
                
                if spC[ihi,ili]=="":###把空的换成"0"
                    spCdata[ihi,ili]=0
                else:
                    spCdata[ihi,ili]=float(spC[ihi,ili])
                ili +=1
            ihi +=1
        CLLXarr.append(spCdata)
    elif signGYL==True:###处理供应量数据
        spJ=np.array(inputarray[1:arrrow,start[0]+1:start[0]+1+datalen])
        hi=spJ.shape[0]
        li=spJ.shape[1]
        spJdata=np.zeros((hi,li))
##        print(hi,li)
        ihi=0
        while ihi<hi:
            ili=0
            while ili< li:                
                if spJ[ihi,ili]=="":###把空的换成"0"
                    spJdata[ihi,ili]=0
                else:
                    spJdata[ihi,ili]=float(spJ[ihi,ili])
                ili +=1
            ihi +=1
        GYLarr=spJdata[0]###只有一行数据
    elif signJL==True:
##        print("处理距离数据")
        spJ=np.array(inputarray[1:arrrow,start[0]+1:start[0]+1+datalen])
        hi=spJ.shape[0]
        li=spJ.shape[1]
        spJdata=np.zeros((hi,li))
##        print(hi,li)
        ihi=0
        while ihi<hi:
            ili=0
            while ili< li:                
                if spJ[ihi,ili]=="":###把空的换成"0"
                    spJdata[ihi,ili]=0
                else:
                    spJdata[ihi,ili]=float(spJ[ihi,ili])
                ili +=1
            ihi +=1
        JLarr=spJdata
##        print("处理完毕")
##    print("记录处理完毕，清除记录")
    inputstr=[]
    haveRecord=False    
#########============处理记录函数结束


#########读取文档开始#######################################
##############读取过程中也保存了读取的数据##################
##############
def readdata():
##    #########===========因为输入变量都是全局变量所以没有变量传入
##    ##==
##    ######申明全局变量，可以在函数外面改变变量的值########
    global inputstr,spots,sptsarr,XQDZarr,WZLXarr,CLLXarr,JLarr
    global sheet1,nrow,haveRecord
    myrow=0
    moHang=False###末行
    sptsarr=[]
    while myrow <nrow:
    ##    print("第",myrow+1,"行")
        onerow = sheet1.row_values(myrow)#获取行内容    
        if set(onerow)=={""} :###空行
            if haveRecord==True:
                chuliJilu()####处理记录
            else:            
    ##            print("跳过该行")
                myrow +=1
                continue
        else :
            ###读取该行并且加入到inputstr中
            inputstr.append(onerow)
    ##        print("把这一行添加到记录中")
            haveRecord=True
            ####下面是：当读取到末行时也要处理inputstr的数据
            if myrow==nrow-1:
                moHang=True
            if moHang==True:
                if haveRecord==True:                
                    chuliJilu()####处理记录
        myrow +=1
   
#########读取文档结束#######################################    



#################数据转换开始######################################
###############从各个分表中得到物资类型总表####################
####把sptsarr数组转换成矩阵
def sptsarrtoarray():
    global sptsarr
    if type(sptsarr)== type([]):
        sptsarr=np.array(sptsarr)        
    else:
        pass

####得到需求点数据矩阵三维矩阵sptsarr
####sptsarr函数和全局变量sptsarr的名称重复,故函数名改为getsptsarr
def getsptsarr():
    readdata()####从表格中读取数据，更新全局变量
     ####数据处理以下##################    
    sptsarrtoarray()    
    sptsarrAscend=sptsarrascend(spots,sptsarr)###数据按顺序排列
    return sptsarrAscend

###sptsarrascend得到按需求点数字从小到大的顺序排列
def sptsarrascend(spots,sptsarr):
    sptsarrAscend=[]
    spotsInt=spotstoint(spots)
    myspotsInt=sorted(spotsInt)
##    print(spotsInt)
    for i in myspotsInt:
        idx=spotsInt.index(i)
##        print(idx)
        sptsarrAscend.append(sptsarr[idx])
    sptsarrAscend=np.array(sptsarrAscend)
    return sptsarrAscend
####spotstoint得到需求点的数字顺序
def spotstoint(spots):
    spotsInt=[]
    for i in spots:
        Obj=re.findall(r'\d{1}',i)###单个需求点
        spotsInt.append(int(Obj[0]))
    return spotsInt

####cllxtypenum从输入数据中得到车辆类型和数量的表
####输入：inputarray矩阵，start列表，datalen整数常量
def cllxtypenum(inputarray,start,datalen):    
    CLLXtype=[]###车辆类型
    CLLXnum=[]###车辆数量
    shuliangpos=list(inputarray[:,start[0]]).index('辆数')
    CLLXtype=inputarray[0,start[0]+1:start[0]+1+datalen]
    CLLXnumtemp=inputarray[shuliangpos,start[0]+1:start[0]+1+datalen]
    CLLXnum=[int(float(i)) for i in CLLXnumtemp]###i:'1.0'
    CLLXtypenum=dict(zip(CLLXtype, CLLXnum))
    return CLLXtypenum
        
def cllxtitle(inputarray,start,datalen):
    CLLXtitl=inputarray[0,start[0]+1:start[0]+1+datalen]    
    return CLLXtitl

    
######传入需求点需求sptsarr三维矩阵
######传出物资类型总表
def wzlx(sptsarr):
    myWZLX=sptsarr[0,:,0:2]
    return myWZLX

def wzlxzl(myWZLX):
    WZLXzl=myWZLX[:,0]
    return WZLXzl

def wzlxtj(myWZLX):
    WZLXtj=myWZLX[:,1]
    return WZLXtj

####wrtthreeandfour函数把三个矩阵数据写到全局变量中
#########################把六个总结数据写到全局变量中,******不能多次调用该函数，只调用一次
def wrtthreeandsix():###获取最终数据
    global SPTSACDarr,XQDHZarr,WZLXHZarr
    global Zxqdnum,Zxqdzl,Zxqdtj,Zclnum,Zclzl,Zcltj
    global Zgylzl,Zgyltj###供应量重量和体积
    global Zfplzl,Zfpltj###分配量重量和体积
    sptsarrAscend=getsptsarr()###获取按顺序排列的数据
    xqdHZ=xqdhz(sptsarrAscend)###获取需求点汇总数据
    wzlxHZ=wzlxhz(sptsarrAscend,GYLarr)###获取物资类型包含供应量汇总数据
    ####把获取的数据写进全局变量中
    for i in sptsarrAscend:
        SPTSACDarr.append(i)
    SPTSACDarr=np.array(SPTSACDarr)   
    for i in xqdHZ:
        XQDHZarr.append(i)
    XQDHZarr=np.array(XQDHZarr)
    for i in wzlxHZ:
        WZLXHZarr.append(i)
    WZLXHZarr=np.array(WZLXHZarr)
    ###总结数据
    Zxqdnum=len(XQDHZarr[:,0])
    Zxqdzl=round(sum(XQDHZarr[:,-2]),4)
    Zxqdtj=round(sum(XQDHZarr[:,-1]),4)
    Zgylzl=round(sum(wzlxHZ[3,:]*wzlxHZ[0,:]),4)#供应量重量
    Zgyltj=round(sum(wzlxHZ[3,:]*wzlxHZ[1,:]),4)#供应量体积
    Zfplzl=round(sum(wzlxHZ[4,:]*wzlxHZ[0,:]),4)#分配量重量
    Zfpltj=round(sum(wzlxHZ[4,:]*wzlxHZ[1,:]),4)#分配量体积
    ###CLLXarr是一个三维矩阵
    lcl=len(CLLXarr)
    li=0
    while li < lcl:
        Zclnum.append(int(sum(CLLXarr[li][2,:])))
        Zclzl.append(round(sum(CLLXarr[li][2,:]*CLLXarr[li][0,:]),4))
        Zcltj.append(round(sum(CLLXarr[li][2,:]*CLLXarr[li][1,:]),4))
        li+=1

#################数据转换结束######################################  
    
#########得到汇总数据开始##########################################
###########输入各个需求点数据三维矩阵sptsarr
#######输出各个需求点含体积和重量的汇总数据
def xqdhz(sptsarr):
    xqdHZ=[]
    myWZLX=wzlx(sptsarr)##物资类型种类
    WZLXzl =wzlxzl(myWZLX)###种类
    WZLXtj=wzlxtj(myWZLX)###体积
    xqdXQZ=xqdxql(sptsarr)
    for i in xqdXQZ:
        i=list(i)##把数组转换成数组
        spotZL=round(sum(i*WZLXzl),4)
        spotTJ=round(sum(i*WZLXtj),4)
        i.append(spotZL)
        i.append(spotTJ)
        xqdHZ.append(i)
    xqdHZ=np.array(xqdHZ)
    return xqdHZ

###########输入各个需求点数据三维矩阵sptsarr
#######输出各个物资类型含总量的汇总数据
def wzlxhz(sptsarr,GYLarr):
    xqdXQZ=xqdxql(sptsarr)
##    print(xqdXQZ)
    myWZLX=wzlx(sptsarr)##物资类型种类
##    print(myWZLX)
    ###https://docs.scipy.org/doc/numpy-1.10.1/reference/generated/numpy.sum.html
    wzlxh=np.sum(xqdXQZ,axis=0)
    wzlxHZ=list(np.transpose(myWZLX))
    wzlxHZ.append(wzlxh)
    wzlxHZ.append(GYLarr)
    ###增加可分配矩阵和分配系数矩阵
    KFParr=[]
    FPXSarr=[]
    le=len(GYLarr)
    li=0
    while li<le:
        KFParr.append(min(wzlxh[li],GYLarr[li]))
        FPXSarr.append(round(min(wzlxh[li],GYLarr[li])/wzlxh[li],4))
        li+=1
    wzlxHZ.append(KFParr)
    wzlxHZ.append(FPXSarr)
    wzlxHZ=np.array(wzlxHZ)
##    print(wzlxHZ)
    return wzlxHZ

####返回不含统计数据的需求点需求量数据
def xqdxql(sptsarr):
    return sptsarr[:,:,2]###需求量
######打印SPTSACDarr,XQDHZarr,WZLXHZarr,CLLXarr,JLarr
def printall():
    print("="*100)
    print("*"*10,"第一步：读取数据（开始）","*"*10)
    print("需求点分表\n",SPTSACDarr)
    print("需求点汇总表\n",XQDHZarr)
    print("物资类型汇总表\n",WZLXHZarr)
    print("供应量总表\n",GYLarr)
    for idx,x in enumerate(CLLXarr):        
        print("车辆类型表",idx+1,"\n",x)
    print("距离表\n",JLarr)
    print("*"*10,"第一步：读取数据（结束）","*"*10)
    print("="*100)
    
    
#########得到汇总数据结束##########################################

###############测试函数开始########################################
def alltest():
    wrtthreeandsix()
    printall()
###############测试函数结束########################################    

    

####运行alltest()更新全局变量，方便被调用，但是只引入一次并保留这些全局变量
####不要多次引入excelchuli模块，这样会多次运行alltest()，不利于程序优化
####更新：需求点分表\需求点汇总表\物资类型汇总表\车辆类型表\距离表\总需求点数等数据
alltest()
###读文件结束


#####第二步
####根据读取的时间得到最短等待时间路线
#################################################################################################################
#################全局变量：开始#######################################
#################全局变量：开始#######################################
######相当于一个登记区域
###受灾点时间矩阵Tij(点0到点8)
Tij=JLarr
###n:共Zclzl(总车辆重量)只蚂蚁
n=Zclnum###数组
###次数序列
ts=[]
###总等待时间序列
zdts=[]
### Taoij：i到j的信息素浓度
### Taoij初始化为1
Taoij=np.ones([9,9])
Taoij[0,:]=0
minzdt=1e15
minRoute=[]
minRouteJH={}###简化路线
### prevbestDTij：到j的当前最短等待时间
prevbestDTij=np.ones([9,9])
prevbestDTij *= 1e15
### 点和等待时间的关系
spotandtime=[]
######根据蚂蚁数量的不同返回相应的最短时间与最短时间对应的时间
TimeandRoute=[]
###重量阈值系数，体积阈值系数，按体积算(True or False)，剩余空间问题(True or False)
Fenpeixishu=[]
###第一次分配矩阵
Diyicifenpeiarr=[]

###体积阈值
TJyuzhi=0
###重量阈值
ZLyuzhi=0
#################全局变量：结束#######################################
#################全局变量：结束#######################################
##############################优先级Pij函数——开始###########################
###优先级Pij函数youxianjiPij：传入:所有点集合V,已经访问过的点集合yifangwen,概率矩阵Pij
#####################返回：经过处理的Pij
def youxianjiPij(V,yifangwen,Pij):
    weifangwen=V-yifangwen###未访问的点
    if len(weifangwen)>0:###确实有未访问的点，则把已访问过的点的Pij赋值为0
        for i in yifangwen:
            Pij[i,:]=0
    return Pij
##############################优先级Pij函数——结束###########################

##############################更新信息素累计函数——开始######################
###v1.1.13引入——更新函数gengxin():更新信息素的累计
############输入：矩阵DTij
############输入并更改：矩阵Taoij，矩阵XSij以及矩阵prevbestDTij
def gengxin(DTij,XSij):###也会改变DTij,XSij
    global Taoij,prevbestDTij

    ####更新prevbestDTij，开始：
    ###获取有效的DTij元素的位置，非零元素
    (row,col)= DTij.nonzero()###https://blog.csdn.net/u011361880/article/details/73611740
    idxi=len(row)
    myii=0
    while myii < idxi:
        ###根据DTij得到的模型，，保留三位小数
        ###提升目标：该行最小的prevbestDTij值
        tishengmubiao=min(prevbestDTij[row[myii],:])
        ###本次最短:该行最小的DTij的值
        bencizuiduan=nonzeromin(DTij[row[myii],:])
        if bencizuiduan< tishengmubiao:
            ###有提升
            ###达到本次最短
            if DTij[row[myii],col[myii]]==bencizuiduan:
                ##提升率
                if tishengmubiao==1e15:
                    tishenglv=0.05
                else:
                    tishenglv=round((tishengmubiao-bencizuiduan)/tishengmubiao,3)###保留三位小数
                ##更新系数矩阵相应系数
                XSij[row[myii],col[myii]]=5*(1+tishenglv)
                ##更新以前最佳等待时间的矩阵相应等待时间
                prevbestDTij[row[myii],col[myii]]=DTij[row[myii],col[myii]]                
                ##检查Taoij对应的值，挪动原来的值到新地方
                bestTao=max(Taoij[row[myii],:])
                if bestTao==1:
                    pass
                else:
                    ###原来的位置信息素为零
                    posT=list(Taoij[row[myii],:]).index(bestTao)
                    Taoij[row[myii],posT]=0
                Taoij[row[myii],col[myii]]=bestTao
                
            else:
                ###未达到本次最短
                ##系数为0
                XSij[row[myii],col[myii]]=0
                ##Taoij为0
                Taoij[row[myii],col[myii]]=0
        elif bencizuiduan == tishengmubiao:
            ###有保持
            ###达到本次最短
            if DTij[row[myii],col[myii]]==bencizuiduan:
                ##系数为0
                XSij[row[myii],col[myii]]=0
                ##更新以前最佳等待时间的矩阵相应等待时间,可能存在多条最短路线
                prevbestDTij[row[myii],col[myii]]=DTij[row[myii],col[myii]]
                ##检查Taoij对应的值，挪动原来的值到新地方
                bestTao=max(Taoij[row[myii],:])
                Taoij[row[myii],col[myii]]=bestTao
            else:
                ###未达到本次最短
                ##系数为0
                XSij[row[myii],col[myii]]=0
                ##Taoij为0
                Taoij[row[myii],col[myii]]=0
        else:
            ###没有提升也没有保持
            ###判断现在的位置是不是提升目标的位置
            ###若是，则系数为0，Taoij不变
            ###否则系数为0，Taoij为0
            pospre=list(prevbestDTij[row[myii],:]).index(tishengmubiao)
            if pospre==col[myii]:
                pass
            else:
                ##Taoij为0
                Taoij[row[myii],col[myii]]=0
            ##系数为0
            XSij[row[myii],col[myii]]=0
                     
        myii +=1
    prevbestDTij=baoliuyigezuixiao(prevbestDTij)##每行保留一个最小值
    ####更新prevbestDTij，结束。
##############################更新信息素累计函数——结束######################
    
##############################矩阵保留最小值函数——开始######################    
###输入一个矩阵A,每行保留一个最小值
def baoliuyigezuixiao(A):
    lenA=len(A)
    ia=0
    while ia<lenA:
        i=A[ia,:]
        best=min(i)
        leni=len(i)
        ii=0
        while ii<leni:
            if A[ia,ii]==best:
                pass
            else:
                A[ia,ii]=1e15
            ii +=1
        ia +=1
    return A
##############################矩阵保留最小值函数——结束######################

##############################非零最小值函数——开始##########################
###某一行的非零最小值
def nonzeromin(hang):
    mymin=1e15
    for i in hang:
        if i!=0 and i<mymin:
            mymin=i
    return mymin
##############################非零最小值函数——结束##########################

##############################简化路线函数——开始############################
def jianhuaroute(routemin):
    jianhuaRoute=[]
    oneline=[]
    for i in routemin:
        if i==0:
            if oneline != [] and set(oneline)!={0}:
                jianhuaRoute.append(tuple(oneline))
            oneline=[]
        else:
            pass
        oneline.append(i)
    if oneline != [] and set(oneline)!={0}:###处理最后一条路线
        jianhuaRoute.append(tuple(oneline))
    
    jianhuaRoute=list(set(jianhuaRoute))
    ###https://www.w3schools.com/python/ref_list_sort.asp
    jianhuaRoute.sort(key=myfunc)
    jianhuaRoutezd={}
    lenr=len(jianhuaRoute)
    idi=1
    while idi<=lenr:
        jianhuaRoutezd[idi]=jianhuaRoute[idi-1]
        idi+=1
    return jianhuaRoutezd
def myfunc(e):
    return e[1]##0之后的第一个点
##############################简化路线函数——结束############################

##############################寻找路线函数——开始############################
def xunzhaoluxian(Tij,N,n,alpha,belta,Padd,rou):
    global ts,zdts,Taoij,minzdt,minRoute,prevbestDTij
    ###迭代开始
    ###t:迭代次数t-th的初始化
    t=0
    while t<N:
        ###V是所有需要访问的点集合， T 是已经访问过的点集合,
        ###Vk是k-th蚂蚁可以访问的点集合,Vno是Vk中不符合载重的点
        V=[1,2,3,4,5,6,7,8]
        T=[]
        ###总等待时间zdt=0
        zdt=0
        ###Routet:t-th次迭代的路线,t-th迭代的目标函数z值
        Routet=[]
        ###DTij:等待时间矩阵
        DTij=np.zeros([9,9])
        ###XSij:系数矩阵，和DTij相关
        XSij=np.zeros([9,9])
        ### Eataij:i到j的期望
        Eataij=np.ones([9,9])
        Eataij*=.1
        Eataij[0,:]=0
        ##print(Eataij)    
        ###k:k-th蚂蚁的初始化
        k=1
        ###存储临时加入的点
        rk=[]
    ##    print("第",t,"次迭代\n")
        while k<=n:
    ##        print("第",k,"只蚂蚁:")
            ###用一个变量标记在一只蚂蚁内，还是出了一只蚂蚁的范围，方便确定下一只蚂蚁站立的位置。
            theAnt=False        
            Vk=set(V)-set(T)
    ##        print("Vk:\n")
    ##        print(Vk)
            
            ###存储一只蚂蚁的累计时间
            tk=0
            if len(Vk)==0:
                break
            else:
                while True:
    ##                print("theAnt,",theAnt)
                    ###获取当前所在位置lastT,对同一只蚂蚁来说会用到.下一只蚂蚁重新回到起点
                    if theAnt==False:
                        Routet.append(0) #一只蚂蚁的开始点
                        ###lastT:上一个访问的元素，意味着当前的位置在lastT这儿
                        lastT=0                    
                    else:
                        lastT=rk[-1]                    
                        ####r0决定蚂蚁是否继续访问下一个点
                        r0=random.random()
                        if r0 > Padd:###退出
    ##                        print("退出")
                            break
                        else:###继续选择下一个点
    ##                        print("继续")
                            #### lastT的点不能再被访问，故相应的Eata概率为0
                            Eataij[lastT,:]=0
                            ###该点以后也不能被访问了
                            T.append(lastT)
    ##                    
    ##                print("当前位置（第x点）：\n")
    ##                print(lastT)
    ##                print(Eataij)
                    ####如果没有可以选择的点，则退出
                    if sum(Eataij[:,lastT])==0:
    ##                    print("退出")
                        break
                    
                    Pij=np.power(Taoij,alpha)*np.power(Eataij,belta)
    ##                print("Pij:\n",Pij)
                    yifangwen=set(rk)
                    Pij=youxianjiPij(set(V),yifangwen,Pij)###经过处理的Pij
    ##                print("Pij:\n",Pij)
                    sumPij=np.sum(Pij,axis=0)#按列求和
    ##                print("sumPij:\n")
    ##                print(sumPij)
                    sumPij=sumPij[lastT]
    ##                print("sumPij:\n")
    ##                print(sumPij)
                    #SPij从归一化的Pij中取出对应的一列
                    ####如果经过处理后没有可以选择的点，则退出
                    if sumPij==0:
    ##                    print("退出")
                        break
                    Pij=Pij/sumPij
                    SPij=Pij[1:,lastT]
    ##                print("SPij\n")
    ##                print(SPij)
                    ###创建一个包含所有点的字典
                    myd=dict(zip(V,SPij))
    ##                print("myd:\n")
    ##                print(myd)
                    d=sorted(myd.items(),key=lambda item:item[1])# 按值排序，得到d是元组
                    dictd=dict(d) #转换成字典，方便后面累加赋值用
    ##                print("dictd:\n")
    ##                print(dictd)
                    ii=0
                    ###while 完成累加和赋值
                    ###
                    while ii<len(d):
                        if ii==0:
                            dictd[d[ii][0]]=dictd[d[ii][0]]
                        elif ii<len(d)-1:
                            dictd[d[ii][0]]=dictd[d[ii-1][0]]+d[ii][1]
                        else:
                            dictd[d[ii][0]]=1
                        ii +=1
                    #print(dictd)
                    ###转换成元组
                    tupd=tuple(zip(dictd,dictd.values()))
    ##                print(tupd)
                    ###r1是从元组中选择一个点
                    r1=random.random()
                   
    ##                print(r1)
                    ii=0
                    while ii<len(tupd):
                        if r1>tupd[ii][1]:                        
                            ii+=1
                        else:
    ##                        print("选第",tupd[ii][0],"个点")
                            rk.append(tupd[ii][0])
                            Routet.append(tupd[ii][0])
                            tk +=Tij[tupd[ii][0],lastT]
                            DTij[tupd[ii][0],lastT]=tk
                            break
                    theAnt=True
    ##            print("rk\n",rk)
    ##        print("下一只蚂蚁：\n")
            k +=1
        ###Tao增加值的矩阵
        addTao=np.zeros(np.shape(Taoij))
        
        ###查看是否有未访问的点
        weifangwen=set(V)-set(Routet)
        if len(weifangwen)>0:
    ##        print("有未访问的点：",weifangwen)
            addTao *=0
    ##        print(addTao)
        else:
    ##        print("无未访问的点")
            ###计算有效的值
            ts.append(t)
            zdt=sum(sum(DTij))
    ##        print("总时间：",zdt)
            zdts.append(zdt)    
    ##        print("路线Routet:")
    ##        print(Routet)
            if zdt < minzdt:
                minzdt=zdt            
                minRoute=[i for i in Routet]
    ##        XSij=xsjuzhen(DTij)
            gengxin(DTij,XSij)
    ##        print("Taoij\n",Taoij)
    ##        print("prevbestDTij\n",prevbestDTij)
    ##        print("XSij\n",XSij)
            
    ##        print(XSij)
            ########### deltaTaoijt:t-th次结束时更新的信息素增量.
            deltaTaoijt=round(1/np.power(zdt/800,3),3) #根据zdt的值得到模型的值，保留三位小数
    ##        print("deltaTaoijt:",deltaTaoijt)
            addTao=deltaTaoijt*XSij
    ##    print("addTao\n",addTao)
    ##    ### Taoij：i到j的信息素浓度
        Taoij = (1-rou)*Taoij+ addTao
    ##    print("Taoij\n",Taoij)
##        if t==N-1:
##            printroute()###打印重要数据
            
        t +=1
    return jianhuaroute(minRoute)
##############################寻找路线函数——结束############################

##############################打印路线——开始################################
def printroute():
    global Taoij,minzdt,minRoute,minRouteJH
##    print("最终的Taoij:\n",Taoij)
    print('最短总等待时间(min)：',minzdt)
##    print("最短时间对应路线：",minRoute)
    print("最短时间对应简化路线(路线序号:(0点,第一个需求点,第二个需求点,...))：\n",minRouteJH)
##############################打印路线——结束################################
##############################画图函数——开始################################
####画总等待时间
def plotzdt(mayinum,ts,zdts):###传入蚂蚁数量,次数序列,等待时间
    global N,alpha,belta,Padd,rou
    ####图像初始化
    #为了图像中显示中文不乱码，设置下面两个参数
    plt.rcParams['font.sans-serif']=['SimHei']
    plt.rcParams['axes.unicode_minus']=False
    ##fig=plt.figure(figsize=(10,6),dpi=700)
    fig=plt.figure(figsize=(10,6))
    ##sub1=fig.add_subplot(211)#上下两张图
    sub1=fig.add_subplot(111)#单独一张图
    titlestr="迭代次数N:"+str(N)+ ",蚂蚁数n:"+str(mayinum)+",alpha:"+str(alpha)+",belta:"+str(belta)+",每只蚂蚁继续访问概率Padd:"+str(Padd)+",信息素挥发率rou:"+str(rou)
    sub1.plot(ts,zdts, lw=2, ls='-',label="目标函数总等待时间zdt")
    sub1.set_title(titlestr)
    sub1.legend()
    #显示图像
    plt.show()
##############################画图函数——结束################################

###
def zuiduanluxian():
    global Tij,N,n,alpha,belta,Padd,rou
    global ts,zdts,Taoij,minzdt,minRoute,minRouteJH,prevbestDTij
    global TimeandRoute,spotandtime
    TimeandRoute=[]
    spotandtime=[]
    print("!"*10,"请注意:观察最短总等待时间，不合适则重新运行程序","!"*10)
    print("="*100)
    print("*"*10,"第二步：寻找路线（开始）","*"*10)
    for num in n:
        ###初始化开始=======
        ###次数序列
        ts=[]
        ###总等待时间序列
        zdts=[]
        ### Taoij：i到j的信息素浓度
        ### Taoij初始化为1
        Taoij=np.ones([9,9])
        Taoij[0,:]=0
        minzdt=1e15
        minRoute=[]
        minRouteJH={}###简化路线
        ### prevbestDTij：到j的当前最短等待时间
        prevbestDTij=np.ones([9,9])
        prevbestDTij *= 1e15
        onespotandtime=[]
        ###初始化结束=========
        minRouteJH=xunzhaoluxian(Tij,N,num,alpha,belta,Padd,rou)###找到最短路线minRoute
        for i in prevbestDTij:
            onespotandtime.append(min(i))
        spotandtime.append(onespotandtime[1:])
        print("-"*50)
        print("车辆数量",num)
        TimeandRoute.append([minzdt,minRouteJH])
        printroute()
        plotzdt(num,ts,zdts)
        print("-"*50)
    print("*"*10,"第二步：寻找路线（结束）","*"*10)
    print("="*100)
    return TimeandRoute



##############################拷贝cheliangzuhe函数——开始############################
##############################拷贝cheliangzuhe函数——开始############################
###这些函数将根据车辆类型三维数组处理成字典的数组,前提需要获取按体积算的信息

def cllxdicttolist(CLLXtypenum):
    chulishuliang=len(CLLXtypenum)###处理数量
##    print(chulishuliang)
    CLLXTNlist=[]##车辆类型数量数组
    cllxtype= list(CLLXtypenum)
    cllxnum=list(CLLXtypenum.values())
    i=0
    while i < chulishuliang:
        j=0
        while j< cllxnum[i]:
            CLLXTNlist.append(cllxtype[i])
            j+=1
        i +=1
    return CLLXTNlist

###传入车辆类型列表,传出车辆组合列表
def cheliangzuhe(list1):    
    list2 = []
    for i in range(1, len(list1)+1):
        iter1 = itertools.combinations(list1, i)
        list2.append(list(iter1))
    list1=list(itertools.chain.from_iterable(list2))
    list2=[]
    for i in list1:
        if not i in list2:
            list2.append(i)
    return list2

##print(list2)

####蚂蚁值列表value_list = [重量值122，体积值比如23.3,车辆组合比如('车型A', '车型B')]
####蚂蚁数组Mayishuzu=[value_list1,...,value_listn]
####每只蚂蚁=蚂蚁序号：蚂蚁值列表
####蚂蚁字典={每只蚂蚁}={蚂蚁序号1：蚂蚁值列表1，...}

####mayishuzu蚂蚁数组
####车辆组合列表list2,
def mayishuzu(list2,CLLXtitl,CLLXarr):
    Mayishuzu=[]###蚂蚁数组
    for i in list2:
        value_list=[]
        zhongliang=0.0
        tiji=0.0
        for j in i:
            posj=list(CLLXtitl).index(j)
            zhongliang +=CLLXarr[0,posj]
            tiji += CLLXarr[1,posj]
        ###一只蚂蚁的数据value_list
        value_list.append(round(zhongliang,3))
        value_list.append(round(tiji,3))
        value_list.append(i)
        Mayishuzu.append(value_list)
    return  Mayishuzu

#####myszpaixu蚂蚁数组排序
###先按重量升序，再按体积升序
def myszpaixu(Mayishuzu,antijisuan):
    if antijisuan==False:###先考虑按重量
        ###https://www.runoob.com/python/python-func-sorted.html
        ###https://blog.csdn.net/y12345678904/article/details/77507552
        return sorted(Mayishuzu, key=lambda x:(x[0],x[1]))  ###先按重量升序，再按体积升序
    else:##先考虑按体积
        return sorted(Mayishuzu, key=lambda x:(x[1],x[0]))  ###先按体积升序，再按重量升序

#####mayizidian对排好序的蚂蚁给一个编号,返回蚂蚁字典
def mayizidian(Mayishuzupaixu):
    lenm=len(Mayishuzupaixu)
    ii=0
    Mayizidian={}
    while ii<lenm:
        Mayizidian[ii]=Mayishuzupaixu[ii]
        ii+=1
    return Mayizidian

def getmayizidian(CLLXtypenum,CLLXtitl,CLLXarr,antijisuan):
    list1=cllxdicttolist(CLLXtypenum)
    list2=cheliangzuhe(list1)
    Mayishuzu=mayishuzu(list2,CLLXtitl,CLLXarr)
##    print(Mayishuzu)
    Mayishuzupaixu=myszpaixu(Mayishuzu,antijisuan)
##    print(Mayishuzupaixu)
    Mayizidian=mayizidian(Mayishuzupaixu)
##    print(Mayizidian)
    return Mayizidian
def mygetmayizidian(CLLXtypenum,CLLXtitl,CLLXarr,Fenpeixishu):
    Mayizidian=[]
    lc=len(CLLXtypenum)
    li=0
    while li<lc:
        Mayizidian.append(getmayizidian(CLLXtypenum[li],CLLXtitl[li],CLLXarr[li],Fenpeixishu[li][2]))
        li+=1
##    print(Mayizidian)
    return Mayizidian
##############################拷贝cheliangzuhe函数——结束############################
##############################拷贝cheliangzuhe函数——结束############################


def abc(Zfplzl,Zfpltj,Zclzl,Zcltj,gongpingxishu):###Zclzl,Zcltj是数组
    ###Fenpeixishu的一个实例 [[0.7124, 0.4233, True, True], [1.032, 0.6212, True, True]]
    ###重量阈值系数，体积阈值系数，按体积算(True or False)，剩余空间问题(True or False)
    
    Fenpeixishu=[]
    le=len(Zclzl)
    li=0
    while li<le:
        thiszlbz=round((Zclzl[li]/Zfplzl)*gongpingxishu[li],4)
        thistjbz=round((Zcltj[li]/Zfpltj)*gongpingxishu[li],4)
        antijisuan=False
        shengyukongjianwenti=False
        if thiszlbz >= 1 and thistjbz >= 1:
            print("非剩余空间问题")            
        else:
            shengyukongjianwenti=True
            if thiszlbz > thistjbz:##重量比值更大
                antijisuan=True##按体积算
            else:
                antijisuan=False###重量比值更小
        Fenpeixishu.append([thiszlbz,thistjbz,antijisuan,shengyukongjianwenti])
        li+=1
    return Fenpeixishu

def diyicikefenpei(XQDHZarr,WZLXHZarr):##第一次可分配
    xishu=WZLXHZarr[-1]
    xuqiu=XQDHZarr[:,:-2]
    danweizl=WZLXHZarr[0,:]
    danweitj=WZLXHZarr[1,:]
    yicifenpei=[]
    for i in xuqiu:
        hang=i*xishu
        yihang=[int(j) for j in hang]
        zl=sum(yihang*danweizl)
        tj=sum(yihang*danweitj)
        yihang.append(zl)
        yihang.append(tj)
        yicifenpei.append(yihang)
    yicifenpei=np.array(yicifenpei)
    return yicifenpei

###获取第一次分配阈值
def diyicifenpeiyuzhi(Diyicifenpeiarr,Fenpeixishu):
    yicifenpeiyuzhi=[]
    for idx,x in enumerate(Fenpeixishu):
        zlxishu=x[0]
        tjxishu=x[1]
        zlxulie=Diyicifenpeiarr[:,-2]
        tjxulie=Diyicifenpeiarr[:,-1]
        zl=zlxulie*zlxishu
        tj=tjxulie*tjxishu
        d=np.transpose([zl,tj])
##        d = np.hstack((zl,tj))
        yicifenpeiyuzhi.append(d)
    yicifenpeiyuzhi=np.array(yicifenpeiyuzhi)
    return yicifenpeiyuzhi
####获取第一次分配的路线阈值
############################################################
def luxianyuzhi(Diyicifenpeiyuzhi,TimeandRoute,Fenpeixishu):
    Diyiciluxianyuzhi=[]
    for idx,x in enumerate(Diyicifenpeiyuzhi):
        if Fenpeixishu[idx][3]==True: ##是剩余空间分配问题           
            ###处理路线阈值
            Diyiciluxianyuzhi.append(chuliluxianyuzhi(TimeandRoute[idx][1],Diyicifenpeiyuzhi[idx]))
        else:
            print("不用处理路线阈值")
    return Diyiciluxianyuzhi                        

###处理路线阈值
def chuliluxianyuzhi(luxian,yuzhi):
    luxianyuzhi=[]###
    lenr=len(luxian)
    ir=1
    while ir <=lenr:
        mayizl=0.
        mayitj=0.
        lenl=len(luxian[ir])
        il=1###除去0点
        mydian=[]
        while il <lenl:
            mydian.append(luxian[ir][il])
            mayizl+=yuzhi[luxian[ir][il]-1,0]##重量阈值
            mayitj+=yuzhi[luxian[ir][il]-1,1]##体积阈值
            il+=1
        thismayi=[round(mayizl,4),round(mayitj,4)]
        thismayiyudian=[mydian,thismayi]
        luxianyuzhi.append(thismayiyudian)
        ir+=1
    return luxianyuzhi
############################################################

####蚂蚁匹配
############################################################
def mayipipei(Mayiyudian,Mayizidian,CLLXtypenum,Fenpeixishu):
    Mayipipei=[]
    ###非常重要深拷贝
    CLLXtypenum=copy.deepcopy(CLLXtypenum)##深拷贝，参考https://www.runoob.com/w3cnote/python-understanding-dict-copy-shallow-or-deep.html
    for idx,x in enumerate(Fenpeixishu):
##        antijisuan=Fenpeixishu[idx][2]###按体积算
        ###剩余车辆    
##        shengyucheliang=CLLXtypenum[idx]
    ##    print(shengyucheliang)
##        print("剩余车辆：/n",CLLXtypenum[idx])
        mayi=xuanmayi(Fenpeixishu[idx][2],Mayiyudian[idx],Mayizidian[idx],CLLXtypenum[idx])
##        print(mayi)
        ###比以前多了后半部分
        shengyu=mayi[1]
        mayi=jiaohuanmayi(mayi[0],Fenpeixishu[idx][2])
        Mayipipei.append([mayi,shengyu])
    return Mayipipei
###交换蚂蚁
def jiaohuanmayi(Mayixuhao,antijisuan):
    mayixuhao=[]
    if antijisuan==True:##按体积倒序排列
        mayixuhao=sorted(Mayixuhao,key=lambda x:-x[1][1][1])
    else:##按重量倒序排列
        mayixuhao=sorted(Mayixuhao,key=lambda x:-x[1][1][0])
    mayityno=sorted([i[0][1] for i in mayixuhao],reverse=True)
    ###修改蚂蚁序号的类型号，即重新匹配蚂蚁
    for idx,x in enumerate(mayixuhao):
        mayixuhao[idx][0][1]=mayityno[idx]
    ###按照点的顺序重新排序
    mayixuhao.sort(key=xxx)
    return mayixuhao
def xxx(x):
    ###元素x的格式：[[1, 1], [[1], [491.238, 3.097]]]
    return x[1][0][0]
###选蚂蚁            
def xuanmayi(antijisuan,Mayiyudian,Mayizidian,shengyucheliang):
    ###蚂蚁序号
    mayixuhao=[]
##    print("mayi",mayi)
    lma=len(Mayiyudian)
    ima=0
    maxuhao=1
    while ima<lma:
        lenzd=len(Mayizidian)
##        print(lenzd)
        zdi=0
        while zdi<lenzd:
            if antijisuan==True:
                if Mayizidian[zdi][1]<Mayiyudian[ima][1][1]:##比体积
                    pass
                else:
                    ####检查剩余的车辆
                    kexuan=jianchacheliangshengyu(shengyucheliang,Mayizidian[zdi])
##                    print("kexuan",kexuan)
                    if kexuan==True:
                        thism=[[maxuhao,zdi],Mayiyudian[ima]]
                        maxuhao+=1
                        mayixuhao.append(thism)
##                        print(Mayiyudian[ima],"被匹配了")
                        shengyucheliang=gengxinshengyucheliang(shengyucheliang,Mayizidian[zdi])
                        break                      

                    else:
                        pass                
            else:
                if Mayizidian[zdi][0]<Mayiyudian[ima][1][0]:##比重量
                    pass
                else:
                    ####检查剩余的车辆
                    kexuan=jianchacheliangshengyu(shengyucheliang,Mayizidian[zdi])
                    if kexuan==True:
                        thism=[[maxuhao,zdi],Mayiyudian[ima]]
                        maxuhao+=1
                        mayixuhao.append(thism)
##                        print(Mayiyudian[ima],"被匹配了")
                        shengyucheliang=gengxinshengyucheliang(shengyucheliang,Mayizidian[zdi])
                        break 
                    else:
                        pass                
            zdi+=1
##        print(Mayiyudian[ima],"找不到匹配")
        ima+=1
##    print("安排蚂蚁后剩余的车辆(车型号：剩余数量)：",shengyucheliang)
    return [mayixuhao,shengyucheliang]                        

###检查车辆剩余
def jianchacheliangshengyu(shengyucheliang,Mayizidianvals):
    ###车辆组合
    ###https://www.runoob.com/python/att-dictionary-copy.html
    jianchashengyucheliang=shengyucheliang.copy()###非常重要
    zuhe=Mayizidianvals[2]
    for i in zuhe:
##        print(i)
        jianchashengyucheliang[i] -=1
    ####检查剩余数量
    shengyushuliang=jianchashengyucheliang.values()
    for i in shengyushuliang:
        if i <0:
            return False
    return True

####更新剩余车辆
def gengxinshengyucheliang(shengyucheliang,Mayizidianvals):
    ###车辆组合
    zuhe=Mayizidianvals[2]
    for i in zuhe:
        shengyucheliang[i] -=1
##    print(shengyucheliang)
    return shengyucheliang
############################################################

###打印出用到的蚂蚁类型号的具体信息
############################################################
def printusedmayi(Mayipipei,Mayizidian):
    Usedmayi=[]
    print('-'*50)
    print("蚂蚁类型信息([","蚂蚁类型号","，蚂蚁重量","，蚂蚁体积","，车型组合","])：")
    for idx,x in enumerate(Mayipipei):
        print("第",idx+1,"次车辆类型安排：")
        mayi=usedmayi(Mayipipei[idx][0],Mayizidian[idx])
        mayi.sort(key=takefirst)
        for m in mayi:
            print(m)
        Usedmayi.append(mayi)
    
    print("*"*10,"第三步：选择蚂蚁组合（结束）","*"*10)
    print("="*100)
    return Usedmayi
def takefirst(elem):
    return elem[0]
        
def usedmayi(Mayixuhao,Mayizidian):
    Usedmayi=[]
    mayitypenolist=[]
    for i in Mayixuhao:
        mayitypeNO=i[0][1]
        if mayitypeNO in mayitypenolist:
            continue
        else:
            mayitypenolist.append(mayitypeNO)
        chexing=",".join(Mayizidian[mayitypeNO][2])##车型连接,https://www.runoob.com/python/att-string-join.html
        oneline=[mayitypeNO,Mayizidian[mayitypeNO][0],Mayizidian[mayitypeNO][1],chexing]
        Usedmayi.append(oneline)
    return Usedmayi
##    print(mayitypenolist)
############################################################

##############################总配送表函数——开始############################
def zxqdps(XQDHZarr,mayixuhao,Diyicifenpeiyuzhi,spotandtime,TimeandRoute):
    ZxqdPSarr=[]
    for idx,x in enumerate(spotandtime):
        ZxqdPSarr.append(zxqdpsarr(XQDHZarr,mayixuhao[idx][0],Diyicifenpeiyuzhi[idx],spotandtime[idx],TimeandRoute[idx][1]))
    return ZxqdPSarr
def zxqdpsarr(XQDHZarr,mayixuhao,Diyicifenpeiyuzhi,spotandtime,minRouteJH):
    ###########重量阈值、体积阈值、配送蚂蚁序号、配送蚂蚁类型号、配送路线序号、配送时间、满意度
    ZxqdPSarr=np.zeros(tuple(np.array(np.shape(XQDHZarr))+np.array((0,7))))
    #####为了得到ZxqdPSarr
    ##从第xx列开始：
    #####重量阈值、体积阈值、配送蚂蚁序号、配送蚂蚁类型号、配送路线序号、配送时间、满意度
    routenumber=minroutenumber(minRouteJH)
##    print(routenumber)
    xqdzlyuzhi=Diyicifenpeiyuzhi[:,0]###需求点重量阈值
    xqdtjyuzhi=Diyicifenpeiyuzhi[:,1]###需求点体积阈值
    lenxqd=len(xqdtjyuzhi)
    il=0
    lxqd=len(XQDHZarr[0,:])
    while il <lenxqd:
        ZxqdPSarr[il,lxqd]=round(xqdzlyuzhi[il],3)###第lxqd+1列：重量阈值
        ZxqdPSarr[il,lxqd+1]=round(xqdtjyuzhi[il],3)###第+2列：体积阈值
        ###第+3列：配送蚂蚁序号和路线序号一样
        
        ###第+4列：配送蚂蚁类型号
        
        ###第+5列：配送路线序号
        ZxqdPSarr[il,lxqd+4]=routenumber[il]
        ###第+6列：配送时间
        ZxqdPSarr[il,lxqd+5]=spotandtime[il]
        il+=1
    ###第+3列：配送蚂蚁序号        
    ###第+4列：配送蚂蚁类型号
    for i in mayixuhao:
        dian=i[1][0]
        for j in dian:
            ###第+3列：配送蚂蚁序号
            ZxqdPSarr[j-1,lxqd+2]=i[0][0]
            ###第+4列：配送蚂蚁类型号
            ZxqdPSarr[j-1,lxqd+3]=i[0][1]
    return ZxqdPSarr
##############################总配送表函数——结束############################

##########################蚂蚁车型分开函数——开始############################
def mayixhchaifen(Mayixuhao,Mayizidian,CLLXarr,CLLXtitl,WZLXHZarr,Usedmayi):
    mayixhcf=[]
    for idx,x in enumerate(Mayixuhao):
        thismayi=mayixuhaochaifen(Mayixuhao[idx][0],Mayizidian[idx],CLLXarr[idx],CLLXtitl[idx],WZLXHZarr,Usedmayi[idx])
        mayixhcf.append(thismayi)
    return mayixhcf
def mayixuhaochaifen(Mayixuhao,Mayizidian,CLLXarr,CLLXtitl,WZLXHZarr,Usedmayi):
    myxh=[]
    for mayixuhao  in Mayixuhao:
        xuhao=mayixuhao[0]
        typ=xuhao[1]
        mayi=Mayizidian[typ]
        mayichezuhe=list(mayi[2])
        zhongliangtiji=[]
        for che in mayichezuhe:
            p=list(CLLXtitl).index(che)
            zhongliangtiji.append([CLLXarr[0,p],CLLXarr[1,p]])
        myxh.append([xuhao,[mayixuhao[1][0],zhongliangtiji]])
    myxh=mayixuhaowuzifenbiao(myxh,WZLXHZarr,Usedmayi,CLLXarr,CLLXtitl)
    return myxh

### 带有车型和物资数量的蚂蚁序号
def mayixuhaowuzifenbiao(Myxuhaofenkai,WZLXHZarr,Usedmayi,CLLXarr,CLLXtitl):
    wuzifenbiao=[0 for i in WZLXHZarr[2,:]]
    for idx,x in enumerate(Myxuhaofenkai):
        mayityno=x[0][1]
        chexingzuhe=[i[3] for i in Usedmayi if i[0]==mayityno]
        chexingzuhe=chexingzuhe[0]
        chexingzuhe=chexingzuhe.split(',')
##        print(chexingzuhe)
        zhongliangtiji={i:[CLLXarr[0,list(CLLXtitl).index(i)],CLLXarr[1,list(CLLXtitl).index(i)]] for i in chexingzuhe}
##        print(zhongliangtiji)
        chexinghewuzi=[[i,zhongliangtiji[i],wuzifenbiao.copy()] for i in chexingzuhe]###.copy很重要
##        print(chexinghewuzi)
        Myxuhaofenkai[idx][1].pop(1)###原来的重量体积不要了
        Myxuhaofenkai[idx][1].extend([chexinghewuzi])

    return Myxuhaofenkai
##########################蚂蚁车型分开函数——结束############################

##############################路线序号函数——开始############################
def minroutenumber(minRouteJH):
    ####返回到每个点的路线序号
    routenumber=[]
    ###https://blog.csdn.net/xijuezhu8128/article/details/88555161
    routekeys = list(minRouteJH.keys())
    routevalues = list(minRouteJH.values())
##    print(routevalues)
    Dian=dian(minRouteJH)
##    print(Dian)
    maxdian=max(Dian)##最大的点
    idx=0
    lenr=maxdian
    while idx < lenr:
        zhedian=idx+1
        for i in routevalues:
            if zhedian in i:
                routenumber.append(routekeys[routevalues.index(i)])
        idx+=1 
    return routenumber
###根据minRouteJH获取有哪些点
def dian(minRouteJH):
    routevalues = list(minRouteJH.values())
    Dian=[]
    for i in routevalues:
        for j in i:
            Dian.append(j)
    Dian=set(Dian)
    return Dian
##############################路线序号函数——结束############################


##################################################################    
'''根据公平系数得到蚂蚁匹配'''
def gpxsdaomayipipei(Zfplzl,Zfpltj,Zclzl,Zcltj,gongpingxishu,CLLXtypenum,CLLXtitl,CLLXarr,XQDHZarr,WZLXHZarr,TimeandRoute):
    Fenpeixishu=abc(Zfplzl,Zfpltj,Zclzl,Zcltj,gongpingxishu)
##    print(Fenpeixishu)
    ###获取蚂蚁字典的数组##需要
    Mayizidian=mygetmayizidian(CLLXtypenum,CLLXtitl,CLLXarr,Fenpeixishu)
    ##print(Mayizidian)
    ###获取第一次分配的矩阵。可以分配的数量
    ##print("第一次可分配数量：")##需要
    Diyicifenpeiarr=diyicikefenpei(XQDHZarr,WZLXHZarr)
##    print(Diyicifenpeiarr)
    ###获取第一次分配阈值。由体积阈值或者重量阈值确定的必须分配的数量
    ##print("第一次分配阈值：")##需要
    Diyicifenpeiyuzhi=diyicifenpeiyuzhi(Diyicifenpeiarr,Fenpeixishu)
    ##print(Diyicifenpeiyuzhi)
    ####获取第一次分配的路线阈值(蚂蚁与点)
    Diyiciluxianyuzhi=luxianyuzhi(Diyicifenpeiyuzhi,TimeandRoute,Fenpeixishu)
##    print(Diyiciluxianyuzhi)
    ###第一次蚂蚁匹配##需要
    Mayipipei=mayipipei(Diyiciluxianyuzhi,Mayizidian,CLLXtypenum,Fenpeixishu)
##    print(Mayipipei)
    return Mayizidian,Diyicifenpeiarr,Diyicifenpeiyuzhi,Mayipipei
##################################################################

###
##################################################################
'''判断公平系数取值合理与否，
若不合理则更新公平系数并重新进行蚂蚁匹配,直到所有公平系数合理为止,
返回最终合理的公平系数 '''
def gongpingxishuheli(Zfplzl,Zfpltj,Zclzl,Zcltj,gongpingxishu,CLLXtypenum,CLLXtitl,CLLXarr,XQDHZarr,WZLXHZarr,TimeandRoute):
    Mayizidian,Diyicifenpeiarr,Diyicifenpeiyuzhi,Mayipipei=[],[],[],[]
    bianhua=0.1*np.ones(len(gongpingxishu))
    diyicipanduan=[True for i in gongpingxishu]
    heli=[False for i in gongpingxishu]###公平系数合理与否
    zengda=[False for i in gongpingxishu]###是否增大（公平系数）
    print("="*100)
    print("*"*10,"第三步：选择蚂蚁组合（开始）","*"*10)
    print('-'*50)
    print("."*10,"正在寻找合适的公平系数","."*10)
    while sum(heli)<len(gongpingxishu):###直到所有公平系数都合理才退出
        print("当前的公平系数是：",gongpingxishu)
        ###先根据公平系数得到蚂蚁匹配
        Mayizidian,Diyicifenpeiarr,Diyicifenpeiyuzhi,Mayipipei=gpxsdaomayipipei(Zfplzl,Zfpltj,Zclzl,Zcltj,gongpingxishu,CLLXtypenum,CLLXtitl,CLLXarr,XQDHZarr,WZLXHZarr,TimeandRoute)
##        print(Mayipipei)
        ###得到公平系数是否合理
        for idx,x in enumerate(Mayipipei):
            print("第",idx+1,"次车辆类型安排：")
            panduangongpingxishuheli(Mayipipei[idx],TimeandRoute[idx],gongpingxishu,bianhua,diyicipanduan,heli,zengda,idx)
    ##        print("第",idx+1,"次车辆类型安排后，剩余车辆情况：")
    ##        print(Mayipipei[idx][1])
        ###不合理则重新安排蚂蚁
    print('')    
    print("找到合理的公平系数是：",gongpingxishu)
    return gongpingxishu,Mayizidian,Diyicifenpeiarr,Diyicifenpeiyuzhi,Mayipipei
'''panduangongpingxishuheli函数传入公平系数，返回公平系数是否合理
其中涉及到更新公平系数的方法，简单来说公平系数如果不够或多了就变化0.1,如果接近合理值时变化幅度会在上一次基础上减半'''
def panduangongpingxishuheli(Mayipipei,TimeandRoute,gongpingxishu,bianhua,diyicipanduan,heli,zengda,idx):###判断公平系数取值合理与否    
    if len(Mayipipei[0])==len(TimeandRoute[1]):
        heli[idx]=True
    else :
        heli[idx]=False
        print("公平系数过大，导致有",len(TimeandRoute[1])-len(Mayipipei[0]),"条路线没有分配到蚂蚁")
        print("需要减小公平系数")
        if diyicipanduan[idx]:
            ###第一次判断，变化幅度为0.1，即变化幅度保持不变
            pass            
        else:###后续判断
            if bianhua[idx]==.1:##上次变化是0.1
                if zengda[idx]:##前一次增大
                    bianhua[idx]*=.5###幅度减半
                else:#前一次减小
                    ###变化幅度为0.1，即变化幅度保持不变
                    pass
            else:##上次变化小于0.1
                bianhua[idx]*=.5###幅度减半
        diyicipanduan[idx]=False##用完了第一次判断
        zengda[idx]=False###减小公平系数
        gongpingxishu[idx]-=bianhua[idx]###更新下一次迭代使用的公平系数
        return heli
    shengyucheliang= Mayipipei[1]
    cheshu=shengyucheliang.values()
    zongcheshu=sum(list(cheshu))
    if zongcheshu>0:
        heli[idx]=False
        print("公平系数过小，导致有车辆未被分配")
        print(shengyucheliang)
        print("需要增大公平系数")
        if diyicipanduan[idx]:
            ###第一次判断，变化幅度为0.1，即变化幅度保持不变
            pass
        else:##后续判断
            if bianhua[idx]==.1:#上次变化是0.1
                if zengda[idx]:##前一次增大
                    pass#变化幅度保持不变
                else:##前一次减小
                    bianhua[idx]*=.5###幅度减半
            else:##上次变化小于0.1
                bianhua[idx]*=.5###幅度减半
        diyicipanduan[idx]=False##用完了第一次判断               
        zengda[idx]=True###增大公平系数
        gongpingxishu[idx]+=bianhua[idx]###更新下一次迭代使用的公平系数
        return heli
    else:
        heli[idx]=True
    print("公平系数合理")
    return heli

##################################################################

####根据有多个车辆方案，对应于多个需求矩阵，这些需求矩阵都是一样的
##################################################################
def duogesptsacdarr(CLLXarr,SPTSACDarr):
    duogespts=[]
    for idx,x in enumerate(CLLXarr):
        duogespts.append(SPTSACDarr.copy())###.copy很重要
    return duogespts
##################################################################

###第一步分配函数
##################################################################
###第一次配送前修改为可分配数量
def xiugaikefenpeishu(DuogeSPTS,Diyicifenpeiarr):
    fenpeishu=Diyicifenpeiarr[:,:-2]
    for idx,x in enumerate(DuogeSPTS):
        for idy,y in enumerate(x):
            y[:,2]=fenpeishu[idy]
    return DuogeSPTS

###第一步分配函数
def diyibufp(SPTSACDarr,Myxuhaofenkai,ZxqdPSarr,WZLXHZarr):###WZLXHZarr已经发生了些许的变化，增加了几行
    sptacd=[]
    for idx,x in enumerate(SPTSACDarr):
        spts=diyibufenpei(SPTSACDarr[idx],Myxuhaofenkai[idx],ZxqdPSarr[idx],WZLXHZarr)
        sptacd.append(spts)
    return sptacd
def diyibufenpei(SPTSACDarr,Myxuhaofenkai,ZxqdPSarr,WZLXHZarr):
    for idx,x in enumerate(Myxuhaofenkai):
        dian=x[1][0]
        mayixuhao=x[0][0]
        mayitypeno=x[0][1]
        chexingwuzi=x[1][1]
        cheshu=len(chexingwuzi)
        iche=0
        wuzishunxu=[]
        for idn,dn in enumerate(dian):
            ####每个需求点的需求
            Dianarr=SPTSACDarr[dn-1]
            ####物资权重列
            wuziquanzhong=Dianarr[:,4]
            ####物资顺序
            wuzishunxu=[(dn,i,x) for i,x in enumerate(wuziquanzhong)]
            ###物资权重从高到底排序
            wuzishunxu=sorted(wuzishunxu,key=lambda x:-x[2])
            ###点的阈值满足
            lwz=len(WZLXHZarr[1,:])
            yixuantiji=sum(ZxqdPSarr[dn-1,:lwz]*WZLXHZarr[1,:])
            yixuanzhongliang=sum(ZxqdPSarr[dn-1,:lwz]*WZLXHZarr[0,:])
            tijiyuzhi=ZxqdPSarr[dn-1,lwz+3]-yixuantiji
            zhongliangyuzhi=ZxqdPSarr[dn-1,lwz+2]-yixuanzhongliang
            if tijiyuzhi<=0 or zhongliangyuzhi<=0:
                manzuyizhi=True
                break###该点的阈值已经满足,安排下一个点
            else:
                manzuyizhi=False
                
            if iche<cheshu:###还有车辆可以分配
                shengyutiji=chexingwuzi[iche][1][1]-sum(chexingwuzi[iche][2]*WZLXHZarr[1,:])
                shengyuzhongliang=chexingwuzi[iche][1][0]-sum(chexingwuzi[iche][2]*WZLXHZarr[0,:])
##                print(shengyuzhongliang,shengyutiji)
            else :###没有车辆可以分配了
                ####下一只蚂蚁
                break
            for idwz,wz in enumerate(wuzishunxu):
                ###一次分配函数
                wuzifenwan=False###物资分完标志
                wuzifenwan,shengyutiji,shengyuzhongliang,iche,manzuyizhi=yicifenpei(shengyutiji,shengyuzhongliang,wz,ZxqdPSarr,SPTSACDarr,iche,cheshu,wuzifenwan,Myxuhaofenkai,idx,chexingwuzi,WZLXHZarr,manzuyizhi)
                ####本次分配有两种情况：
                ####1.分配达到阈值
                if manzuyizhi:
                    break###退出该点的物资分配，即安排下一个点
                else:
                ####2.没有达到阈值
                    ####2.1.物资分配完了，记录分配到最后的车的序号、剩余空间和剩余体积，以便下一个物资的分配。
                    ####2.2.物资没有分配完，且所有的车辆的空间用完。
                    if wuzifenwan==True:##情况2.1.
                        continue##下一种物资分配
                    else:###情况2.2.
                        break###退出该点的物资分配，即安排下一只蚂蚁
    ###更新配送表中的体积和重量
    lwz=len(WZLXHZarr[0,:])#物资个数
    for idx,x in enumerate(ZxqdPSarr):
        ZxqdPSarr[idx,lwz+1]=sum(ZxqdPSarr[idx,:lwz]*WZLXHZarr[1,:])##体积
        ZxqdPSarr[idx,lwz]=sum(ZxqdPSarr[idx,:lwz]*WZLXHZarr[0,:])##重量
    return SPTSACDarr   

def yicifenpei(shengyutiji,shengyuzhongliang,wz,ZxqdPSarr,SPTSACDarr,iche,cheshu,wuzifenwan,Myxuhaofenkai,idx,chexingwuzi,WZLXHZarr,manzuyizhi):
    dianpos=wz[0]-1
    wuzipos=wz[1]
    Dianarr=SPTSACDarr[dianpos]
    ###车辆限制，在体积和重量的限制下的可选物资数量,向下取整
    shu1=int(shengyutiji/Dianarr[wuzipos,1])##体积数
    shu2=int(shengyuzhongliang/Dianarr[wuzipos,0])##重量数
    chezai_shuliang=min(shu1,shu2)
    ##点的阈值，在体积阈值和重量阈值的限制下的可选物资数量，向上取整
    lwz=len(WZLXHZarr[0,:])#物资个数
    yixuantiji=sum(ZxqdPSarr[dianpos,:lwz]*WZLXHZarr[1,:])
    yixuanzhongliang=sum(ZxqdPSarr[dianpos,:lwz]*WZLXHZarr[0,:])
    tijiyuzhi=ZxqdPSarr[dianpos,lwz+3]-yixuantiji
    zhongliangyuzhi=ZxqdPSarr[dianpos,lwz+2]-yixuanzhongliang
    tijishu=math.ceil(tijiyuzhi/Dianarr[wuzipos,1])
    zhongliangshu=math.ceil(zhongliangyuzhi/Dianarr[wuzipos,0])
    dian_shuliang=min(tijishu,zhongliangshu)
    ###点的物资待配送数量
    daisong_shuliang=int(Dianarr[wuzipos,2]-Dianarr[wuzipos,3])
    if chezai_shuliang >=daisong_shuliang:###车辆能够满足待配送物资
        if daisong_shuliang<dian_shuliang:###未达到阈值
            ###该物资完全配送
            ZxqdPSarr[dianpos,wuzipos]+=daisong_shuliang
            ###可以更改SPTSACDarr的第四列####
            SPTSACDarr[dianpos,wuzipos,3]+=daisong_shuliang
            ####更改蚂蚁序列中相应的蚂蚁车型物资数量
            Myxuhaofenkai[idx][1][1][iche][2][wuzipos]+=daisong_shuliang
            ###更新每个需求点剩余的体积和重量阈值
            shengyutiji-=Dianarr[wuzipos,1]*daisong_shuliang
            shengyuzhongliang-=Dianarr[wuzipos,0]*daisong_shuliang
            wuzifenwan=True
            manzuyizhi=False
        else:##满足或超过阈值
            ###该物资完全配送
            ZxqdPSarr[dianpos,wuzipos]+=dian_shuliang
            ###可以更改SPTSACDarr的第四列####
            SPTSACDarr[dianpos,wuzipos,3]+=dian_shuliang
            ####更改蚂蚁序列中相应的蚂蚁车型物资数量
            Myxuhaofenkai[idx][1][1][iche][2][wuzipos]+=dian_shuliang
            ###更新每个需求点剩余的体积和重量阈值
            shengyutiji-=Dianarr[wuzipos,1]*dian_shuliang
            shengyuzhongliang-=Dianarr[wuzipos,0]*dian_shuliang
            if daisong_shuliang==dian_shuliang:
                wuzifenwan=True
            else:
                wuzifenwan=False
            manzuyizhi=True
            
    else:###车辆不能够满足待配送物资
        ###该物资部分配送
        if chezai_shuliang<dian_shuliang:###未达到阈值
            ZxqdPSarr[dianpos,wuzipos]+=chezai_shuliang
            ###可以更改SPTSACDarr的第四列####
            SPTSACDarr[dianpos,wuzipos,3]+=chezai_shuliang
            ####更改蚂蚁序列中相应的蚂蚁车型物资数量
            Myxuhaofenkai[idx][1][1][iche][2][wuzipos]+=chezai_shuliang
            ###更新每个需求点剩余的体积和重量阈值
            shengyutiji-=Dianarr[wuzipos,1]*chezai_shuliang
            shengyuzhongliang-=Dianarr[wuzipos,0]*chezai_shuliang
            wuzifenwan=False
            manzuyizhi=False
        else:##满足或超过阈值
            ZxqdPSarr[dianpos,wuzipos]+=dian_shuliang
            ###可以更改SPTSACDarr的第四列####
            SPTSACDarr[dianpos,wuzipos,3]+=dian_shuliang
            ####更改蚂蚁序列中相应的蚂蚁车型物资数量
            Myxuhaofenkai[idx][1][1][iche][2][wuzipos]+=dian_shuliang
            ###更新每个需求点剩余的体积和重量阈值
            shengyutiji-=Dianarr[wuzipos,1]*dian_shuliang
            shengyuzhongliang-=Dianarr[wuzipos,0]*dian_shuliang
            wuzifenwan=False
            manzuyizhi=True
            
    if manzuyizhi:###满足阈值，该点不需要继续分配了
        return wuzifenwan,shengyutiji,shengyuzhongliang,iche,manzuyizhi
    else:###没有满足阈值，需要继续分配        
        ###物资没有分配完,继续分配物资
        if wuzifenwan==False:
            ####车辆空间不够了       
            ###检查还有没有可以分配的车
            iche+=1
            if iche<cheshu:
                ###上次剩余的空间不考虑了
                shengyutiji=chexingwuzi[iche][1][1]-sum(chexingwuzi[iche][2]*WZLXHZarr[1,:])
                shengyuzhongliang=chexingwuzi[iche][1][0]-sum(chexingwuzi[iche][2]*WZLXHZarr[0,:])
                ###进行下一次的一次分配
                return yicifenpei(shengyutiji,shengyuzhongliang,wz,ZxqdPSarr,SPTSACDarr,iche,cheshu,wuzifenwan,Myxuhaofenkai,idx,chexingwuzi,WZLXHZarr,manzuyizhi)
            else:
                ###车辆用完了
                return wuzifenwan,shengyutiji,shengyuzhongliang,iche,manzuyizhi
            ####车辆空间充足,进行分配
            
        else:###物资分配完了
            return wuzifenwan,shengyutiji,shengyuzhongliang,iche,manzuyizhi
###第一次分配后剩余的可分配物资个数
def diyicifenpeishengyu(ZxqdPSarr,WZLXHZarr):
    fenpei=[]
    kefenpei=WZLXHZarr[4,:]
    for idx,x in enumerate(ZxqdPSarr):
        thisx=sum(ZxqdPSarr[idx][:,:len(kefenpei)])
        thx=[int(i-j) for i,j in zip(kefenpei,thisx)]
        fenpei.append(thx)
    return fenpei                            
##################################################################

###第二步分配，利用蚂蚁剩余空间
##################################################################
###第二次分配前恢复需求矩阵中的各点需求量
def huifuxuqiuliang(DuogeSPTS,SPTSACDarr):
    for idx,x in enumerate(DuogeSPTS):
        for idy,y in enumerate(x):
            y[:,2]=SPTSACDarr[idy][:,2]
    return DuogeSPTS
######第二步分配
def dierbufp(SPTSACDarr,Myxuhaofenkai,ZxqdPSarr,WZLXHZarr,Fenpeishengyu):
    sptsa=[]
    for idx,x in enumerate(SPTSACDarr):
        spts=dierbufenpei(SPTSACDarr[idx],Myxuhaofenkai[idx],ZxqdPSarr[idx],WZLXHZarr,Fenpeishengyu[idx])
        sptsa.append(spts)
    return sptsa    
def dierbufenpei(SPTSACDarr,Myxuhaofenkai,ZxqdPSarr,WZLXHZarr,Fenpeishengyu):
    for idx,x in enumerate(Myxuhaofenkai):
        dian=x[1][0]
        mayixuhao=x[0][0]
        mayitypeno=x[0][1]
        zhongliangtijizuhe=x[1][1]
        cheshu=len(zhongliangtijizuhe)
        iche=0
        wuzishunxu=[]
        for idn,dn in enumerate(dian):
            ####每个需求点的需求
            Dianarr=SPTSACDarr[dn-1]
            ####物资权重列
            wuziquanzhong=Dianarr[:,4]
            ####物资顺序
            wuzishunxu.extend([(dn,i,x) for i,x in enumerate(wuziquanzhong)])
        ###每只蚂蚁的需求物资排好序        
        ###物资权重从高到底排序
        wuzishunxu=sorted(wuzishunxu,key=lambda x:-x[2])
            
        if iche<cheshu:###还有车辆可以分配
            shengyutiji=zhongliangtijizuhe[iche][1][1]-sum(zhongliangtijizuhe[iche][2]*WZLXHZarr[1,:])
            shengyuzhongliang=zhongliangtijizuhe[iche][1][0]-sum(zhongliangtijizuhe[iche][2]*WZLXHZarr[0,:])
        else :###没有车辆可以分配了
            ####下一只蚂蚁
            break
        ###开始权重从高到底开始分配####
        for idwz,wz in enumerate(wuzishunxu):
            ###一次分配函数
            wuzifenwan=False###物资分完标志
            wuzifenwan,shengyutiji,shengyuzhongliang,iche=yicifenpei_2(shengyutiji,shengyuzhongliang,wz,ZxqdPSarr,SPTSACDarr,iche,cheshu,wuzifenwan,Myxuhaofenkai,idx,zhongliangtijizuhe,WZLXHZarr,Fenpeishengyu)
            ####本次分配有两种情况：
            ####1.物资分配完了，记录分配到最后的车的序号，以及剩余体积和剩余重量，以便下一个物资的分配
            ####2.物资没有分配完，且所有的车辆的空间用完。
            if wuzifenwan==True:##情况1
                continue##下一种物资分配
            else:###情况2
                break###退出该点的物资分配，即安排下一只蚂蚁
    ###更新配送表中的体积和重量
    lwz=len(WZLXHZarr[1,:])
    for idx,x in enumerate(ZxqdPSarr):
        ZxqdPSarr[idx,lwz+1]=sum(ZxqdPSarr[idx,:lwz]*WZLXHZarr[1,:])
        ZxqdPSarr[idx,lwz]=sum(ZxqdPSarr[idx,:lwz]*WZLXHZarr[0,:])
    return SPTSACDarr
def yicifenpei_2(shengyutiji,shengyuzhongliang,wz,ZxqdPSarr,SPTSACDarr,iche,cheshu,wuzifenwan,Myxuhaofenkai,idx,chexingwuzi,WZLXHZarr,Fenpeishengyu):
    dianpos=wz[0]-1
    wuzipos=wz[1]
    Dianarr=SPTSACDarr[dianpos]
    ###车辆限制，在体积和重量的限制下的可选物资数量,向下取整
    shu1=int(shengyutiji/Dianarr[wuzipos,1])##体积数
    shu2=int(shengyuzhongliang/Dianarr[wuzipos,0])##重量数
    chezai_shuliang=min(shu1,shu2)
    ###点的物资待配送数量，受到可分配物资的制约
    daisong_shuliang=int(min(Dianarr[wuzipos,2]-Dianarr[wuzipos,3],Fenpeishengyu[wuzipos]))
    if chezai_shuliang >=daisong_shuliang:###车辆能够满足待配送物资
        ###该物资完全配送
        ZxqdPSarr[dianpos,wuzipos]+=daisong_shuliang
        ###更改SPTSACDarr的第四列####
        SPTSACDarr[dianpos,wuzipos,3]+=daisong_shuliang
        ####更改蚂蚁序列中相应的蚂蚁车型物资数量
        Myxuhaofenkai[idx][1][1][iche][2][wuzipos]+=daisong_shuliang
        ####更新可分配剩余量
        Fenpeishengyu[wuzipos]-=daisong_shuliang
        ###更新每个需求点剩余的体积和重量阈值
        shengyutiji-=Dianarr[wuzipos,1]*daisong_shuliang
        shengyuzhongliang-=Dianarr[wuzipos,0]*daisong_shuliang
        wuzifenwan=True
    else:###车辆不能够满足待配送物资
        ZxqdPSarr[dianpos,wuzipos]+=chezai_shuliang
        ###更改SPTSACDarr的第四列####
        SPTSACDarr[dianpos,wuzipos,3]+=chezai_shuliang
        ####更改蚂蚁序列中相应的蚂蚁车型物资数量
        Myxuhaofenkai[idx][1][1][iche][2][wuzipos]+=chezai_shuliang
        ####更新可分配剩余量
        Fenpeishengyu[wuzipos]-=chezai_shuliang
        ###更新每个需求点剩余的体积和重量阈值
        shengyutiji-=Dianarr[wuzipos,1]*chezai_shuliang
        shengyuzhongliang-=Dianarr[wuzipos,0]*chezai_shuliang
        wuzifenwan=False
    ###若物资没有分配完,则继续分配物资
    if wuzifenwan==False:
        ####车辆空间不够了       
        ###检查还有没有可以分配的车
        iche+=1
        if iche<cheshu:
            ###上次剩余的空间不考虑了
            shengyutiji=chexingwuzi[iche][1][1]-sum(chexingwuzi[iche][2]*WZLXHZarr[1,:])
            shengyuzhongliang=chexingwuzi[iche][1][0]-sum(chexingwuzi[iche][2]*WZLXHZarr[0,:])
            ###进行下一次的一次分配
            return yicifenpei_2(shengyutiji,shengyuzhongliang,wz,ZxqdPSarr,SPTSACDarr,iche,cheshu,wuzifenwan,Myxuhaofenkai,idx,chexingwuzi,WZLXHZarr,Fenpeishengyu)
        else:
            ###车辆用完了
            return wuzifenwan,shengyutiji,shengyuzhongliang,iche
        ####车辆空间充足,进行分配
        
    else:###物资分配完了
        return wuzifenwan,shengyutiji,shengyuzhongliang,iche
##################################################################
    
###报告初次分配的信息
##################################################################
def diyicifenpeizongjie(ZxqdPSarr,Myxuhaofenkai,Fenpeishengyu,WZLXHZarr):
    wuzishu=len(Fenpeishengyu[0])
    print("="*100)
    print("*"*10,"第四步：初次分配（开始）","*"*10)
    for idx,x in enumerate(ZxqdPSarr):
        print('-'*50)
        print("按方案",idx+1,"分配,分配矩阵为：")
        print(" 下表标题：1物资1配送数量","，2物资2配送数量","，3物资配送3数量","，4物资4配送数量",'，5配送重量',"，6配送体积",'，7重量阈值','，8体积阈值','，9蚂蚁序号','，10蚂蚁类型号','，11路线序号','，12等待时间(min)','，13满意度(暂未定义)')
        print(ZxqdPSarr[idx])
        print("")
        print("按方案",idx+1,"分配,装车情况为：")
##        print("下表元素格式：[[蚂蚁序号,蚂蚁类型号][[需求点号],[[车型1,[车型1的重量，车型1的体积],[车型1装载物资1的个数...]][车型2...]...]]]")
##        print(Myxuhaofenkai[idx])
        prntmayishengyu(Myxuhaofenkai[idx],WZLXHZarr)
        print("")
        print("按方案",idx+1,"分配后,物资1到物资",wuzishu,"的剩余可分配量")
        print(Fenpeishengyu[idx])
    print("*"*10,"第四步：初次分配（结束）","*"*10)
    print("="*100)
    
###汇总统计和报告最终分配的信息
##################################################################
def pshuizong(ZxqdPSarr,WZLXHZarr):
    pshz=[]
    for idx,x in enumerate(ZxqdPSarr):
        pshz.append(list(peisonghuizong(ZxqdPSarr[idx],WZLXHZarr)))
    return pshz
def peisonghuizong(ZxqdPSarr,WZLXHZarr):
    wuzishu=len(WZLXHZarr[0,:])
    wzzltj=ZxqdPSarr[:,:wuzishu+2]
    return sum(wzzltj)
def printzongjie(DuogeSPTS,ZxqdPSarr,Myxuhaofenkai,Zxqdnum,Zclnum,Zxqdzl,Zxqdtj,Zgylzl,Zgyltj,Zfplzl,Zfpltj,Zclzl,Zcltj,Fenpeishengyu,PShz,WZLXHZarr):
    print("="*100)
    print("*"*10,"第五步：最终分配（开始）","*"*10)
    print('-'*50)
    print('*'*5,'汇总信息','*'*5)
    print("总需求点数：",Zxqdnum)
    for idx,x in enumerate(Zclnum):
        print("方案",idx+1,"的总车辆数：",Zclnum[idx])
        
    print("总需求点重量：",Zxqdzl,"总需求点体积：",Zxqdtj)
    print("总供应量重量：",Zgylzl,"总供应量体积：",Zgyltj)
    print("总分配量重量：",Zfplzl,"总分配量体积：",Zfpltj)
    wuzishu=len(WZLXHZarr[0,:])###物资数量
    lz=len(Zclnum)
    li=0
    while li<lz:
        print("方案",li+1,"：总车辆重量：",Zclzl[li],"总车辆体积：",Zcltj[li])
        print("方案",li+1,"：总配送重量：",round(PShz[li][wuzishu],4),"总配送体积：",round(PShz[li][wuzishu+1],4))
        li+=1
    print("需求点物资1到物资",wuzishu,"的需求量")
    xqllist=[int(i) for i in WZLXHZarr[2,:]]
    print(xqllist)
    print("配送中心物资1到物资",wuzishu,"的供应量")
    gyllist=[int(i) for i in WZLXHZarr[3,:]]
    print(gyllist)
    print("进行分配时物资1到物资",wuzishu,"的可分配量")
    kfplist=[int(i) for i in WZLXHZarr[4,:]]
    print(kfplist)
    for idx,x in enumerate(Fenpeishengyu):
        print("")
        print("按方案",idx+1,"分配,物资1到物资",wuzishu,"的实际分配量")
        sjfpl=[int(i) for i in PShz[idx][:wuzishu]]
        print(' ',sjfpl)
        print("按方案",idx+1,"分配后,物资1到物资",wuzishu,"的剩余可分配量")
        syfpl=[int(i) for i in Fenpeishengyu[idx]]
        print(syfpl)
    print('-'*50)
    print('-'*50)
    print('*'*5,'具体方案','*'*5)
    for idx,x in enumerate(ZxqdPSarr):
        
        print("按方案",idx+1,"分配,分配矩阵为：")
        print(" 下表标题：1物资1配送数量","，2物资2配送数量","，3物资配送3数量","，4物资4配送数量",'，5配送重量',"，6配送体积",'，7重量阈值','，8体积阈值','，9蚂蚁序号','，10蚂蚁类型号','，11路线序号','，12等待时间(min)','，13满意度(暂未定义)')
        print(ZxqdPSarr[idx])
        print("")
        print("按方案",idx+1,"分配,装车情况为：")
##        print("下表元素格式：[[蚂蚁序号,蚂蚁类型号][[需求点号],[[车型1,[车型1的重量，车型1的体积],[车型1装载物资1的个数...]][车型2...]...]]]")
##        print(Myxuhaofenkai[idx])
        prntmayishengyu(Myxuhaofenkai[idx],WZLXHZarr)
        print("")
        print("按方案",idx+1,"分配,需求点的需求矩阵为：")
        print(" 下表的标题为：单位重量，单位体积，需求个数，配送个数，物资权重")
        print(DuogeSPTS[idx])
        print("")
    print("*"*10,"第五步：最终分配（结束）","*"*10)
    print("="*100)

def prntmayishengyu(Myxuhaofenkai,WZLXHZarr):
    danweizl=WZLXHZarr[0,:]
    danweitj=WZLXHZarr[1,:]
    for idx,x in enumerate(Myxuhaofenkai):
        cheliang=x[1][1]
        print(" 蚂蚁序号或路线号",x[0][0],"按顺序经过",x[1][0],'点，车辆的空间剩余情况和物资装载情况：')
        for idy,y in enumerate(cheliang):
            chex=y[0]
            chexzltj=y[1]
            shiyongzl=sum(y[2]*danweizl)
            shiyongtj=sum(y[2]*danweitj)
            shengyuzl=round(chexzltj[0]-shiyongzl)
            shengyutj=round(chexzltj[1]-shiyongtj,2)
##            if shengyuzl>0 and shengyutj>0:                
            print('  ',chex,"剩余重量：",shengyuzl,"剩余体积：",shengyutj,'运送物资个数：',y[2])
            
        
        
##################################################################        
####
####
#########手动设置参数区域#######################
###路线相关参数----------------------------
###最大迭代次数N
N=8000
###alpha:信息素浓度Taoij的权值
alpha=6
###belta：期望Eataij的权值
belta=1
###每只蚂蚁继续访问的概率
Padd=0.2
###信息素挥发率
rou=0.
###----------------------------
###蚂蚁匹配相关参数----------------------------
####设置公平系数
gongpingxishu=0.7*np.ones(len(Zclzl),dtype=float)##初始化的公平系数都为0.7
#########手动设置参数区域#########################

###获取最短时间和路线的组合数组
TimeandRoute=zuiduanluxian()
##print(TimeandRoute)
##print(spotandtime)
###获取分配系数的数组

##gongpingxishu需要根据mayipipei打印出来的结果进行调整。
##################################################################
gongpingxishu,Mayizidian,Diyicifenpeiarr,Diyicifenpeiyuzhi,Mayipipei=gongpingxishuheli(Zfplzl,Zfpltj,Zclzl,Zcltj,gongpingxishu,CLLXtypenum,CLLXtitl,CLLXarr,XQDHZarr,WZLXHZarr,TimeandRoute)
####

###把用到的蚂蚁的组合类型打印出来
Usedmayi=printusedmayi(Mayipipei,Mayizidian)
####获取总需求点分配数组
ZxqdPSarr=zxqdps(XQDHZarr,Mayipipei,Diyicifenpeiyuzhi,spotandtime,TimeandRoute)
Myxuhaofenkai=mayixhchaifen(Mayipipei,Mayizidian,CLLXarr,CLLXtitl,WZLXHZarr,Usedmayi)
##print(Myxuhaofenkai)
##
#####第一次分配
####根据有多个车辆方案，对应于多个需求矩阵，这些需求矩阵都是一样的
DuogeSPTS=duogesptsacdarr(CLLXarr,SPTSACDarr)
##print(DuogeSPTS)
###修改可分配数量
DuogeSPTS=xiugaikefenpeishu(DuogeSPTS,Diyicifenpeiarr)
##print(DuogeSPTS)
##print(SPTSACDarr)
DuogeSPTS=diyibufp(DuogeSPTS,Myxuhaofenkai,ZxqdPSarr,WZLXHZarr)
##print(DuogeSPTS)
###第一次分配剩余
Fenpeishengyu=diyicifenpeishengyu(ZxqdPSarr,WZLXHZarr)
###报告初次分配的信息
diyicifenpeizongjie(ZxqdPSarr,Myxuhaofenkai,Fenpeishengyu,WZLXHZarr)
DuogeSPTS=huifuxuqiuliang(DuogeSPTS,SPTSACDarr)
##print(DuogeSPTS)        
###第二次分配
DuogeSPTS=dierbufp(DuogeSPTS,Myxuhaofenkai,ZxqdPSarr,WZLXHZarr,Fenpeishengyu)
##print(DuogeSPTS)
##print(ZxqdPSarr)
##print(Fenpeishengyu)
PShz=pshuizong(ZxqdPSarr,WZLXHZarr)##
##print(PShz)
###汇总统计和报告最终分配的信息
printzongjie(DuogeSPTS,ZxqdPSarr,Myxuhaofenkai,Zxqdnum,Zclnum,Zxqdzl,Zxqdtj,Zgylzl,Zgyltj,Zfplzl,Zfpltj,Zclzl,Zcltj,Fenpeishengyu,PShz,WZLXHZarr)

