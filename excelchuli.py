####读写Excel文件
####当前版本：v2.1.5
######能从分表中直接读取获得总表数据，存到全局变量中，并且打印出来
#####引用：https://zhuanlan.zhihu.com/p/38492442
#####改进1：增加记录车辆车型和个数的矩阵CLLXtypenum
#####问题：第二次用不了
import xlrd#读Excle
import numpy as np
import re###正则表达式


##取消科学计数法https://blog.csdn.net/Andrew_jdw/article/details/82350041
np.set_printoptions(suppress=True)


#################全局变量：开始#######################################
#################全局变量：开始#######################################
####引用：https://www.runoob.com/python/python-reg-expressions.html
####引用：https://blog.csdn.net/u010412858/article/details/83062200
wb = xlrd.open_workbook("更新表格信息.xlsx")#打开文件
##print(wb.sheet_names())#获取所有表格名字
sheet1 = wb.sheet_by_index(1)#通过索引获取表格

############################################
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
Zclnum=0###总车辆个数
Zclzl=0###总车辆重量
Zcltj=0###总车辆体积
############################################
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
############################################

############################################
###处理距离矩阵要用到的变量###
patternJL = re.compile(r'距离')
###距离表
JLarr=[]##############################################重要
############################################

############################################
#####总表开始
#####需求点汇总
XQDHZarr=[]###########################################重要
#####物资类型汇总
WZLXHZarr=[]##########################################重要
#####总表结束
############################################



#####暂时不会用到！
############################################
###处理需求点总表要用到的变量###
patternZ = re.compile(r' *需求点$')
XQDZarr=[]###需求点总表
############################################

#####暂时不会用到！
############################################
###处理应急物资类型要用到的变量###
patternWZLX = re.compile(r' *应急物资类型')
WZLXarr=[]###应急物资类型表
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
    
    #########==========输入变量有：全局变量有：pattern,patternZ,patternWZLX,patternCLLX,patternJL
    #########======================全局变量需要改变的有：haveRecord,inputstr,spots,sptsarr,XQDZarr,WZLXarr,CLLXarr,JLarr
    #########===========因为输入变量都是全局变量所以没有变量传入
    ##==
    ######申明全局变量，可以在函数外面改变变量的值########
    global pattern,patternZ,patternWZLX,patternCLLX,patternJL
    global haveRecord,inputstr,spots,sptsarr,XQDZarr,WZLXarr,CLLXarr,JLarr
    global CLLXtypenum,CLLXtitl
        
    ###处理类型标志初始化#####
    signXQD=False ###处理需求点1到需求点n的标志
    signXQDZ=False ###处理需求点总表的标志
    signWZLX=False ###处理应急物资类型的标志
    signCLLX=False ###处理车辆类型的标志
    signJL=False ###处理距离矩阵的标志
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
                objXQDZ=patternZ.findall(i)###需求点汇总                        
                objWZLX = patternWZLX.findall(i)###物资类型
                objCLLX=patternCLLX.findall(i)###车辆类型
                objJL=patternJL.findall(i)###距离
                if Obj !=[]:
##                    print(Obj)
                    signXQD=True
                    spots.append(Obj[0])###将需求点加入到顺序表中
                    start.append(idx)###将开始位置记录下来
                if objXQDZ !=[]:
##                    print(objXQDZ)
                    signXQDZ=True
                    start.append(idx)###将开始位置记录下来
                if objWZLX !=[]:
##                    print(objWZLX)
                    signWZLX=True
                    start.append(idx)###将开始位置记录下来
                if objCLLX !=[]:
##                    print(objCLLX)
                    signCLLX=True

                    start.append(idx)###将开始位置记录下来
                if objJL !=[]:
##                    print(objJL)
                    signJL=True

                    start.append(idx)###将开始位置记录下来
                
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
                        spotdata[ihi,ili]=0
                    else:
                        spotdata[ihi,ili]=float(spot[ihi,ili])
                    ili +=1
                ihi +=1
            sptsarr.append(spotdata)###把需求点数据矩阵保存好
            sti +=1
##        print("处理完毕")
    elif signXQDZ==True:
##        print("处理需求点总表数据")
        XQDZarr=inputarray[1:arrrow,start[0]+1:start[0]+1+datalen]###需求总表数据
##        print("处理完毕")
    elif signWZLX==True:
##        print("处理物资类型数据")
        WZLXarr=inputarray[1:arrrow,start[0]+1:start[0]+1+datalen]###物资类型数据
##        print("处理完毕")
    elif signCLLX==True:
##        test1forcllx(inputarray,start,datalen)
        CLLXtypenum=cllxtypenum(inputarray,start,datalen)
        CLLXtitl=cllxtitle(inputarray,start,datalen)
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
        CLLXarr=spCdata
##        print("处理完毕")
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
##    print(sptsarr)
##    print(type(sptsarr))

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
        
def test1forcllx(inputarray,start,datalen):
    print(inputarray)
    print(start[0])
    print(datalen)

###
##def cllxtitle(inputarray,start,datalen,CLLXtitl):
##    CLLXtitltemp=inputarray[0,start[0]+1:start[0]+1+datalen]
##    for i in CLLXtitltemp:
##        CLLXtitl.append(i)
####    print(CLLXtitl)
##    return CLLXtitl

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
def wrtthreeandsix():
    global SPTSACDarr,XQDHZarr,WZLXHZarr
    global Zxqdnum,Zxqdzl,Zxqdtj,Zclnum,Zclzl,Zcltj
    
    sptsarrAscend=getsptsarr()###获取按顺序排列的数据
    xqdHZ=xqdhz(sptsarrAscend)###获取需求点汇总数据
    wzlxHZ=wzlxhz(sptsarrAscend)###获取物资类型汇总数据
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
    Zxqdnum=len(XQDHZarr[:,4])
    Zxqdzl=round(sum(XQDHZarr[:,4]),4)
    Zxqdtj=round(sum(XQDHZarr[:,5]),4)  
    Zclnum=int(sum(CLLXarr[2,:]))
    Zclzl=round(sum(CLLXarr[2,:]*CLLXarr[0,:]),4)
    Zcltj=round(sum(CLLXarr[2,:]*CLLXarr[1,:]),4)
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
##    print(xqdXQZ)
##    print(myWZLX)
##    print(WZLXzl)
##    print(WZLXtj)
    for i in xqdXQZ:
        i=list(i)##把数组转换成数组
##        print(i)
##        print(i*WZLXzl)
        spotZL=round(sum(i*WZLXzl),4)
##        print(spotZL)
        spotTJ=round(sum(i*WZLXtj),4)
##        print(spotTJ)
        i.append(spotZL)
        i.append(spotTJ)
##        print(i)
        xqdHZ.append(i)
    xqdHZ=np.array(xqdHZ)
    return xqdHZ

###########输入各个需求点数据三维矩阵sptsarr
#######输出各个物资类型含总量的汇总数据
def wzlxhz(sptsarr):
    xqdXQZ=xqdxql(sptsarr)
##    print(xqdXQZ)
    myWZLX=wzlx(sptsarr)##物资类型种类
##    print(myWZLX)
    ###https://docs.scipy.org/doc/numpy-1.10.1/reference/generated/numpy.sum.html
    wzlxh=np.sum(xqdXQZ,axis=0)
    wzlxHZ=list(np.transpose(myWZLX))
    wzlxHZ.append(wzlxh)
    wzlxHZ=np.array(wzlxHZ)
##    print(wzlxHZ)
    return wzlxHZ

####返回不含统计数据的需求点需求量数据
def xqdxql(sptsarr):
    return sptsarr[:,:,2]###需求量

##########一个函数获取所有需要的信息
def getall(sptsarrAcd,xqdHZ,wzlxHZ):
    sptsarrAcd=getsptsarr()
    xqdHZ=xqdhz(sptsarrAcd)
    wzlxHZ=wzlxhz(sptsarrAcd)

######打印SPTSACDarr,XQDHZarr,WZLXHZarr,CLLXarr,JLarr
def printall():
    print("需求点分表\n",SPTSACDarr)
    print("需求点汇总表\n",XQDHZarr)
    print("物资类型汇总表\n",WZLXHZarr)
    print("车辆类型表\n",CLLXarr)
    print("距离表\n",JLarr)
    print("总需求点数：",Zxqdnum,"总需求点重量：",Zxqdzl,"总需求点体积：",Zxqdtj)
    print("总车辆数：",Zclnum,"总车辆重量：",Zclzl,"总车辆体积：",Zcltj)
#########得到汇总数据结束##########################################

###############测试函数开始########################################
def xqdhztest():
    sptsarr=getsptsarr()
    xqdHZ=xqdhz(sptsarr)
    print(xqdHZ)

def wzlxhztest():
    sptsarr=getsptsarr()
    wzlxHZ=wzlxhz(sptsarr)
    print(wzlxHZ)
def alltest():
    wrtthreeandsix()
##    printall()
###############测试函数结束########################################    

    

####运行alltest()更新全局变量，方便被调用，但是只引入一次并保留这些全局变量
####不要多次引入excelchuli模块，这样会多次运行alltest()，不利于程序优化
####更新：需求点分表\需求点汇总表\物资类型汇总表\车辆类型表\距离表\总需求点数等数据
##alltest()

