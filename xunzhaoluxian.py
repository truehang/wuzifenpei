#################################################################################################################
#####xunzhaoluxian模块说明：#####################################################################################
#####模块名称：xunzhaoluxian 寻找路线
#####主要特点：寻找最短等待时间路线---收敛有时会比较慢，但最迟会在不超过1万次时收敛，一般在2000次以内收敛
#####继承关系：代码来自lunwen1.1.14.py,对过程进行了函数化，方便被调用
#####主要子区域：1、引入模块区域；2、函数区域；3、全局变量区域；4、函数测试区域
#####版本号：v1.11
#####更新过程：
#############更新1：简化重复的路线
#############更新2：引入配送环节
#########################体积比值和重量比值，分析哪一个成为限制因素。判断是否是满意度问题
#########################路线序号
#############更新3：引入蚂蚁字典
#############更新4：根据###按体积算antijisuan=False?来确定对蚂蚁字典的排序方式
#############更新5：公平系数rougongpingxishurou=0.6###0.9、0.8、0.7不行.
#############更新6：增加蚂蚁序号和类型号,并更新分配总表
#############更新7：第一步分配，每个需求点满足阈值，权重值从高到底满足。
#############更新8：获取蚂蚁剩余空间，引入需求点待运物资数组存储供剩余空间运输的剩余物资
#############更新9：第二步分配，利用蚂蚁剩余空间
#############更新10：把用到的蚂蚁的组合类型打印出来
#############更新11：改进全局变量的赋值方式，安排一个函数来赋值
#####问题：
#############问题1：会发生最短总等待时间增加了。最短总等待时间变成： 720.0
#############问题2：在jupyter上运行只能运行第一次
#####解决问题：
#############问题1：目前没有完全解决，只能重新运行一次观察时间是否是705
#####
#################################################################################################################

#################################################################################################################
###################################引入模块区域##################################################################
import random
import numpy as np
import matplotlib.pyplot as plt #绘制图像要用到的pyplot
import excelchuli as exch###引入自己编写的excelchuli.py
import itertools
import math
##取消科学计数法https://blog.csdn.net/Andrew_jdw/article/details/82350041
np.set_printoptions(suppress=True)
###################################引入模块区域##################################################################
#################################################################################################################


#################################################################################################################
###################################函数区域######################################################################

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
##############################寻找路线函数——结束############################

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

##############################路线序号函数——开始############################
def minroutenumber(minRouteJH):
    ####返回到每个点的路线序号
    routenumber=[]
    ###https://blog.csdn.net/xijuezhu8128/article/details/88555161
    routekeys = list(minRouteJH.keys())
    routevalues = list(minRouteJH.values())
    idx=0
    lenr=8
    while idx < lenr:
        zhedian=idx+1
        for i in routevalues:
            if zhedian in i:
                routenumber.append(routekeys[routevalues.index(i)])
        idx+=1 
    return routenumber
##############################路线序号函数——结束############################

###蚂蚁阈值
def mayiyuzhi(ZxqdPSarr,minRouteJH):
    Mayiyudian=[]###蚂蚁与点，格式：[[点],[蚂蚁的重量，体积]]
    lenr=len(minRouteJH)
    ir=1
##    mayi=[]
    
    while ir <=lenr:
        mayizl=0.
        mayitj=0.
        lenl=len(minRouteJH[ir])
        il=1##除去0点
        ###蚂蚁与点
        mayiyudian=[]
        while il <lenl:
            mayiyudian.append(minRouteJH[ir][il])
            mayizl+=ZxqdPSarr[minRouteJH[ir][il]-1,6]###第7列：重量阈值
            mayitj+=ZxqdPSarr[minRouteJH[ir][il]-1,7]###第8列：体积阈值
            il+=1
        thismayi=[round(mayizl,3),round(mayitj,3)]
        thismayiyudian=[mayiyudian,thismayi]
        Mayiyudian.append(thismayiyudian)
        ir +=1
    return Mayiyudian

##############################总配送表函数——开始############################
def zxqdpsarr():
    #####为了得到ZxqdPSarr
    ##从第7列开始：
    #####重量阈值、体积阈值、配送蚂蚁序号、配送蚂蚁类型号、配送路线序号、配送时间、满意度
    global ZxqdPSarr,XQDHZarr
    global TJyuzhi,ZLyuzhi
    global prevbestDTij,minRouteJH
    
    routenumber=minroutenumber(minRouteJH)
    xqdzlyuzhi=XQDHZarr[:,-2]*ZLyuzhi###需求点重量阈值
    xqdtjyuzhi=XQDHZarr[:,-1]*TJyuzhi###需求点体积阈值
    lenxqd=len(xqdtjyuzhi)
    il=0
    while il <lenxqd:
        ZxqdPSarr[il,6]=round(xqdzlyuzhi[il],3)###第7列：重量阈值
        ZxqdPSarr[il,7]=round(xqdtjyuzhi[il],3)###第8列：体积阈值
        ###第9列：配送蚂蚁序号
        
        ###第10列：配送蚂蚁类型号
        
        ###第11列：配送路线序号
        ZxqdPSarr[il,10]=routenumber[il]
        ###第12列：配送时间
        ZxqdPSarr[il,11]=min(prevbestDTij[il+1])
        il+=1
##############################总配送表函数——开始############################
        
##############################更新问题分类——开始############################
def updateclass():
    global zlbz,tjbz,manyiduwenti,antijisuan
    if zlbz>1 and tjbz>1:
        manyiduwenti=False
    else:
        manyiduwenti=True
    if tjbz<=zlbz:
        antijisuan=True
    else:
        antijisuan=False
    
##############################更新问题分类——结束############################

##############################拷贝cheliangzuhe函数——开始############################
##############################拷贝cheliangzuhe函数——开始############################
###车辆类型数量字典转成列表
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
            tiji += exch.CLLXarr[1,posj]
        ###一只蚂蚁的数据value_list
        value_list.append(round(zhongliang,3))
        value_list.append(round(tiji,3))
        value_list.append(i)
        Mayishuzu.append(value_list)
    return  Mayishuzu

#####myszpaixu蚂蚁数组排序
###先按重量升序，再按体积升序
def myszpaixu(Mayishuzu):
    global antijisuan###考虑按体积算这个因素
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

###一个函数得到蚂蚁字典
def getmayizidian():
    list1=cllxdicttolist(CLLXtypenum)
    list2=cheliangzuhe(list1)
    Mayishuzu=mayishuzu(list2,CLLXtitl,CLLXarr)
##    print(Mayishuzu)
    Mayishuzupaixu=myszpaixu(Mayishuzu)
##    print(Mayishuzupaixu)
    Mayizidian=mayizidian(Mayishuzupaixu)
##    print(Mayizidian)
    return Mayizidian
##############################拷贝cheliangzuhe函数——结束############################
##############################拷贝cheliangzuhe函数——结束############################

####蚂蚁匹配
def mayipipei(Mayiyudian,Mayizidian,CLLXtypenum):
    global antijisuan###按体积算
    ###剩余车辆    
    shengyucheliang=CLLXtypenum
##    print(shengyucheliang)
    mayixuhao=xuanmayi(antijisuan,Mayiyudian,Mayizidian,shengyucheliang)
    return mayixuhao
    

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
                if Mayizidian[zdi][0]<mayi[ima][1][0]:##比重量
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
    print("安排蚂蚁后剩余的车辆(车型号：剩余数量)：",shengyucheliang)
    return mayixuhao
                        

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

###更新总分配表的蚂蚁序号
def updatezfpbmyxh(mayixuhao):
    ###mayixuhao的每个元素形似：[[3, 4], [[3, 5], [715.949, 4.53]]]
    #####################解释：[[蚂蚁序号,蚂蚁类型号],[[点1,点2],[蚂蚁重量，蚂蚁体积]]]
    global ZxqdPSarr
    ###第9列：配送蚂蚁序号        
    ###第10列：配送蚂蚁类型号
    for i in mayixuhao:
        dian=i[1][0]
        for j in dian:
            ###第9列：配送蚂蚁序号
            ZxqdPSarr[j-1,8]=i[0][0]
            ###第10列：配送蚂蚁类型号
            ZxqdPSarr[j-1,9]=i[0][1]
        
##########################第一步分配函数——开始##############################
###第一步分配函数
def diyibufenpei():
    global SPTSACDarr,ZxqdPSarr
    lenSP=len(SPTSACDarr)
    dth=0
    while dth<lenSP:
        ###每个需求点的需求
        Dianarr=SPTSACDarr[dth]
        ####物资权重列
        wuziquanzhong=Dianarr[:,4]
        ####物资顺序
        wuzishunxu=[]
        ii=0
        for i in wuziquanzhong:
            wuzishunxu.append((ii,i))
            ii+=1
        ###物资权重从高到底排序
        wuzishunxu=sorted(wuzishunxu,key=lambda x:x[1],reverse=True)
##        print(wuzishunxu)
        ###每个需求点剩余的体积和重量阈值
        shengyutiji=ZxqdPSarr[dth,7]
        shengyuzhongliang=ZxqdPSarr[dth,6]
        ###开始权重从高到底开始分配#####################################################
        ###hi第几种物资，权重从高到低
        hi=0
        lenh=len(wuzishunxu)###共有几种物资        
        while hi<lenh:
            ####当前物资位置
            wuzipos=wuzishunxu[hi][0]
            ###当下物资可选体积数量和重量数量
            shu1=math.ceil(shengyutiji/Dianarr[wuzipos,1])##体积数
            shu2=math.ceil(shengyuzhongliang/Dianarr[wuzipos,0])##重量数            
            if min(shu1,shu2)>=Dianarr[wuzipos,2]:
                ###该物资完全配送
                ZxqdPSarr[dth,wuzipos]=Dianarr[wuzipos,2]
                ###可以更改SPTSACDarr的第四列，全局变量####
                SPTSACDarr[dth,wuzipos,3]=Dianarr[wuzipos,2]
                ###更新每个需求点剩余的体积和重量阈值
                shengyutiji-=Dianarr[wuzipos,1]*Dianarr[wuzipos,2]
                shengyuzhongliang-=Dianarr[wuzipos,0]*Dianarr[wuzipos,2]
##                print("分配完",wuzishunxu[hi],"后，剩余体积：",shengyutiji,"剩余重量:",shengyuzhongliang)
            else:
                ###该物资部分配送
                ZxqdPSarr[dth,wuzipos]=min(shu1,shu2)
                ###可以更改SPTSACDarr的第四列，全局变量####
                SPTSACDarr[dth,wuzipos,3]=min(shu1,shu2)
                ###更新每个需求点剩余的体积和重量阈值
                shengyutiji-=Dianarr[wuzipos,1]*min(shu1,shu2)
                shengyuzhongliang-=Dianarr[wuzipos,0]*min(shu1,shu2)
##                print("分配到",wuzishunxu[hi],"剩余体积：",shengyutiji,"剩余重量:",shengyuzhongliang)
            ####检查剩余重量和剩余体积
            if shengyutiji<=0 or shengyuzhongliang<=0:##等号很重要
                ###阈值已经分配完毕，退出                
##                print("阈值已经分配完毕，退出")
                break
            hi +=1
        ###每个需求点按阈值分配完物资，开始汇总物资################################
        fenpeishuliang=ZxqdPSarr[dth,:lenh]
        feipeitiji=sum(Dianarr[:,1]*fenpeishuliang)##分配体积
        feipeizhongliang=sum(Dianarr[:,0]*fenpeishuliang)##分配重量
        ZxqdPSarr[dth,5]=feipeitiji
        ZxqdPSarr[dth,4]=feipeizhongliang
        dth +=1



##########################第一步分配函数——结束##############################

        
######################获取蚂蚁剩余空间函数——开始############################
def mayishengyukongjian():
    global ZxqdPSarr,Mayixuhao,Mayizidian,Mayishengyukongjian
    Mayishengyukongjian=[]
    for i in Mayixuhao:
        ###在一只蚂蚁内
        ###蚂蚁类型号
        mayitypeNO=i[0][1]
##        print("mayitypeNO",mayitypeNO)
        yixuantj=0.##已选体积
        yixuanzl=0.##已选重量
        xuandian=i[1][0]
        for j in xuandian:
            yixuantj+=ZxqdPSarr[j-1,5]##将需求点分派体积加上
            yixuanzl+=ZxqdPSarr[j-1,4]##将需求点分派重量加上
##        print(Mayizidian[mayitypeNO][1],Mayizidian[mayitypeNO][0])
##        print(yixuantj,yixuanzl)        
        shengyutj=round(Mayizidian[mayitypeNO][1]-yixuantj,3)##剩余体积
        shengyuzl=round(Mayizidian[mayitypeNO][0]-yixuanzl,3)##剩余重量
##        print(shengyutj,shengyuzl)
        if shengyutj<0 or shengyuzl<0:
            print("蚂蚁安排出错！请重新安排蚂蚁。")
            break
        zhezhimayi=[i[0],[xuandian,[shengyuzl,shengyutj]]]
        Mayishengyukongjian.append(zhezhimayi)
        
######################获取蚂蚁剩余空间函数——结束############################

######################获取蚂蚁剩余空间函数——结束############################

        
########################需求点待运物资函数——开始############################
def xqddaiyunwuzi():
    global SPTSACDarr,Xqddaiyunwuzi
    Xqddaiyunwuzi={}###点序号：物资排序
    lenSP=len(SPTSACDarr)
    dth=0
    while dth<lenSP:
        ###每个需求点的需求
        Dianarr=SPTSACDarr[dth]
        lenD=len(Dianarr)
        ii=0
        yigexqddaiyun=[]###一个需求点待运物资
        while ii<lenD:
            yihang=Dianarr[ii]##一行
            if yihang[3]<yihang[2]:
                yigexqddaiyun.append([dth+1,[ii,yihang[4]]])
            ii+=1        
        ###物资权重从高到底排序
        yigexqddaiyun=sorted(yigexqddaiyun,key=lambda x:x[1],reverse=True)
        Xqddaiyunwuzi[dth+1]=yigexqddaiyun
        dth+=1
########################需求点待运物资函数——结束############################
        
############################第二步分配函数——开始############################
###第二步分配，利用蚂蚁剩余空间,测试 完成
def dierbufenpei():
    global ZxqdPSarr,SPTSACDarr,Xqddaiyunwuzi,Mayishengyukongjian
    ###对每只蚂蚁进行判断
    for i in Mayishengyukongjian:
        shengxiazl=i[1][1][0]##剩余重量
        shengxiatj=i[1][1][1]##剩余体积
        ###待处理点集
        dianlist=i[1][0]
        ###待处理点物资
        chulidianwuzi=[]
        for j in dianlist:
            wuzilist=Xqddaiyunwuzi[j]
            for wuzi in wuzilist:
                chulidianwuzi.append(wuzi)
        chulidianwuzi=sorted(chulidianwuzi,key=lambda x:(-x[1][1],x[0]))
##        print(chulidianwuzi)
        ###结合待处理点需求表来处理物资
        for ich in chulidianwuzi:
            dth=ich[0]-1
            hangpos=ich[1][0]
            zhedian=SPTSACDarr[dth]
            zlshu=math.ceil(shengxiazl/zhedian[hangpos,0])##按重量算的数量
            tjshu=math.ceil(shengxiatj/zhedian[hangpos,1])##按体积算的数量
            shengyu=zhedian[hangpos,2]-zhedian[hangpos,3]
            if min(zlshu,tjshu)>=shengyu:
                ZxqdPSarr[dth,hangpos]+=shengyu
                SPTSACDarr[dth,hangpos,3]+=shengyu
                shengxiazl-=zhedian[hangpos,0]*shengyu
                shengxiatj-=zhedian[hangpos,1]*shengyu
            else:
                ZxqdPSarr[dth,hangpos]+=min(zlshu,tjshu)
                SPTSACDarr[dth,hangpos,3]+= min(zlshu,tjshu)
                shengxiazl-=zhedian[hangpos,0]*min(zlshu,tjshu)
                shengxiatj-=zhedian[hangpos,1]*min(zlshu,tjshu)
            if shengxiatj<=0 or shengxiazl<=0:
##                print("该需求点处理完毕")
                break
        if shengxiazl>0 and shengxiatj>0:
            print("注意：蚂蚁序号为",i[0][0],"的蚂蚁未充分利用完剩余空间。")
    ####
    wuzidanweizl=SPTSACDarr[0][:,0]
    wuzidanweitj=SPTSACDarr[0][:,1]
    lenwuzi=len(wuzidanweitj)##
    lendian=len(ZxqdPSarr)##
    il=0
    while il < lendian:
        ZxqdPSarr[il,4]=sum(wuzidanweizl*ZxqdPSarr[il,:lenwuzi])
        ZxqdPSarr[il,5]=sum(wuzidanweitj*ZxqdPSarr[il,:lenwuzi])
        il+=1
            
############################第二步分配函数——结束############################

##########################打印用到蚂蚁信息——开始############################
###打印出用到的蚂蚁类型号的具体信息
def usedmayi():
    global Mayixuhao,Mayizidian,Usedmayi
    Usedmayi=[]
    mayitypenolist=[]
    for i in Mayixuhao:
        mayitypeNO=i[0][1]
        if mayitypeNO in mayitypenolist:
            continue
        else:
            mayitypenolist.append(mayitypeNO)
        chexing=",".join(Mayizidian[mayitypeNO][2])##车型连接,https://www.runoob.com/python/att-string-join.html
        oneline=[str(mayitypeNO),str(Mayizidian[mayitypeNO][0]),str(Mayizidian[mayitypeNO][1]),chexing]
        Usedmayi.append(oneline)   
##    print(mayitypenolist)
    

##########################打印用到蚂蚁信息——结束############################        
    
######################更新全局变量函数——开始################################
def updatequanjubianliang():
    global Tij,SPTSACDarr,XQDHZarr,WZLXHZarr,CLLXarr,CLLXtitl,CLLXtypenum,JLarr
    global Zxqdnum,Zxqdzl,Zxqdtj,Zclnum,Zclzl,Zcltj
    ###从excel文件中读取
    ####运行alltest()更新excelchuli全局变量,再把全局变量读取过来保存
    exch.alltest()
    ###这次用到的变量
    ###受灾点时间矩阵Tij(点0到点8),矩阵元素之间不能用空格符隔开,
    ###用逗号隔开，时间单位min
    Tij=exch.JLarr

    ####后面会用的变量，从excel中读取的全局变量

    ###需求点分表
    SPTSACDarr=exch.SPTSACDarr
    ###需求点汇总表
    XQDHZarr=exch.XQDHZarr
    ###物资类型汇总表
    WZLXHZarr=exch.WZLXHZarr
    ###车辆类型表
    CLLXarr=exch.CLLXarr
    ###车辆类型标题
    CLLXtitl=exch.CLLXtitl
    ###车辆类型标题
    CLLXtypenum=exch.CLLXtypenum
    ###距离表
    JLarr=exch.JLarr
    ###总需求点数
    Zxqdnum=exch.Zxqdnum
    ###总需求点重量
    Zxqdzl=exch.Zxqdzl
    ###总需求点体积
    Zxqdtj=exch.Zxqdtj
    ###总车辆数
    Zclnum=exch.Zclnum
    ###总车辆重量
    Zclzl=exch.Zclzl
    ###总车辆体积
    Zcltj=exch.Zcltj
######################更新全局变量函数——结束################################

##############################示例函数——开始################################
##############################示例函数——结束################################

##############################画图函数——开始################################
####画总等待时间
def plotzdt():
    global N,n,alpha,belta,Padd,rou,ts,zdts
    ####图像初始化
    #为了图像中显示中文不乱码，设置下面两个参数
    plt.rcParams['font.sans-serif']=['SimHei']
    plt.rcParams['axes.unicode_minus']=False
    ##fig=plt.figure(figsize=(10,6),dpi=700)
    fig=plt.figure(figsize=(10,6))
    ##sub1=fig.add_subplot(211)#上下两张图
    sub1=fig.add_subplot(111)#单独一张图
    titlestr="迭代次数N:"+str(N)+ ",蚂蚁数n:"+str(n)+",alpha:"+str(alpha)+",belta:"+str(belta)+",每只蚂蚁继续访问概率Padd:"+str(Padd)+",信息素挥发率rou:"+str(rou)
    sub1.plot(ts,zdts, lw=2, ls='-',label="目标函数总等待时间zdt")
    sub1.set_title(titlestr)
    sub1.legend()
    #显示图像
    plt.show()
##############################画图函数——结束################################

#################保存变量(可以补充)

##############################打印函数——开始################################
#####打印excel中的数据
def printall():
    print("需求点分表\n",SPTSACDarr)
    print("需求点汇总表\n",XQDHZarr)
    print("物资类型汇总表\n",WZLXHZarr)
    print("车辆类型表\n",CLLXarr)
    print("距离表\n",JLarr)
    print("总需求点数：",Zxqdnum,"总需求点重量：",Zxqdzl,"总需求点体积：",Zxqdtj)
    print("总车辆数：",Zclnum,"总车辆重量：",Zclzl,"总车辆体积：",Zcltj)


def printroute():
    global Taoij,minzdt,minRoute,minRouteJH
##    print("最终的Taoij:\n",Taoij)
    print('最短总等待时间(min)：',minzdt)
##    print("最短时间对应路线：",minRoute)
##    minRouteJH=jianhuaroute(minRoute)###更新minRouteJH
    print("最短时间对应简化路线(路线序号:(0点,第一个需求点,第二个需求点,...))：\n",minRouteJH)
##############################打印函数——结束################################



###################################函数区域######################################################################
#################################################################################################################



#################################################################################################################
###################################全局变量区域##################################################################

#########条件输入区域#############################################################
#####全局变量移到updatequanjubianliang中
updatequanjubianliang()
####后面会用的，从excel中读取的全局变量

#########条件输入区域#############################################################

#########手动设置参数区域#########################################################
###最大迭代次数N
N=8000
###n:共Zclzl(总车辆重量)只蚂蚁
n=Zclnum
###alpha:信息素浓度Taoij的权值
alpha=6
###belta：期望Eataij的权值
belta=1
###每只蚂蚁继续访问的概率
Padd=0.2
rou=0.
#########手动设置参数区域#########################################################

#########其余全局变量区域#########################################################
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

##########物资分配用到的全局变量——开始####################
###重量比值
zlbz=round(Zclzl/Zxqdzl,3)
###体积比值
tjbz=round(Zcltj/Zxqdtj,3)
###按体积算
antijisuan=False
###满意度问题
manyiduwenti=False
###公平系数rou
gongpingxishurou=0.6###0.9、0.8、0.7不行.
###体积阈值
TJyuzhi=tjbz*gongpingxishurou
###重量阈值
ZLyuzhi=zlbz*gongpingxishurou
###需求点配送总表,在XQDHZarr的基础上增加7列：
###########重量阈值、体积阈值、配送蚂蚁序号、配送蚂蚁类型号、配送路线序号、配送时间、满意度
ZxqdPSarr=np.zeros(tuple(np.array(np.shape(XQDHZarr))+np.array((0,7))))
##########物资分配用到的全局变量——结束####################

       

##########蚂蚁字典用到的全局变量——开始####################
Mayizidian={}##蚂蚁字典
Mayixuhao=[]##蚂蚁序号
##Mayiyudian=[]###蚂蚁与点，格式：[[点],[蚂蚁的重量，体积]]
##########蚂蚁字典用到的全局变量——结束####################

##########分配时用到的全局变量——开始######################
###蚂蚁剩余空间
Mayishengyukongjian=[]
###需求点待运物资
Xqddaiyunwuzi={}
##########分配时用到的全局变量——结束######################

####打印格式调整
Usedmayi=[]


#########其余全局变量区域#########################################################

###################################全局变量区域##################################################################
#################################################################################################################

#################################################################################################################
###################################函数测试区域##################################################################

xunzhaoluxian(Tij,N,n,alpha,belta,Padd,rou)###找到最短路线minRoute
minRouteJH=jianhuaroute(minRoute)###更新minRouteJH
printroute()
print("\n")
#####如果最短路线不对的话重新运行上面三条语句
plotzdt()
##print(Tij)
##printall()

###########--------------------------------------------------------------------------############################
##############分配
zxqdpsarr()
##print(sum(ZxqdPSarr[:,9]))
##print(ZxqdPSarr)
Mayiyudian=mayiyuzhi(ZxqdPSarr,minRouteJH)
##print(Mayiyudian)
###########--------------------------------------------------------------------------############################
##############蚂蚁字典
updateclass()###更新问题分类变量
##print("zlbz",zlbz,"tjbz",tjbz,"antijisuan",antijisuan,"manyiduwenti",manyiduwenti)
Mayizidian=getmayizidian()
##print(Mayizidian)
##print(CLLXtypenum)
Mayixuhao=mayipipei(Mayiyudian,Mayizidian,CLLXtypenum)
print("蚂蚁序号信息（[[蚂蚁序号，蚂蚁类型号],[[需求点号],[需求点总重量阈值，需求点总体积阈值]]]）：")
print(Mayixuhao)
###把用到的蚂蚁的组合类型打印出来
usedmayi()
print("蚂蚁类型信息([","蚂蚁类型号","，蚂蚁重量","，蚂蚁体积","，车型组合","])：")
print(Usedmayi)
updatezfpbmyxh(Mayixuhao)##更新总分配表的蚂蚁号
##print(ZxqdPSarr)
###########--------------------------------------------------------------------------############################
##############第一步分配
###作为对比，先打印出需求点的体积重量
##print(XQDHZarr)
###根据阈值进行第一步分配
diyibufenpei()
##print(ZxqdPSarr)
##print(SPTSACDarr)
###########--------------------------------------------------------------------------############################
###得到蚂蚁剩余空间
mayishengyukongjian()
##print(Mayishengyukongjian)
####引入需求点待运物资数组存储供剩余空间运输的剩余物资，增加xqddaiyunwuzi()函数
xqddaiyunwuzi()
##print(Xqddaiyunwuzi)
###########--------------------------------------------------------------------------############################
###第二步分配，利用蚂蚁剩余空间
dierbufenpei()
print("\n")
print("需求点分配方案开始：",'='*70)
print("下表标题：1物资1配送数量","，2物资2配送数量","，3物资配送3数量","，4物资4配送数量",'，5配送重量',"，6配送体积",'，7重量阈值','，8体积阈值','，9蚂蚁序号','，10蚂蚁类型号','，11路线序号','，12等待时间(min)','，13满意度(暂未定义)')
print(ZxqdPSarr)
print('下表标题：需求点物资单位重量','，需求点物资单位体积','，需求点物资数量','，需求点物资配送量','，物资权重')
print(SPTSACDarr)
print("需求点分配方案完",'='*70)
###########--------------------------------------------------------------------------############################
######打印信息
#####把用到的蚂蚁的组合类型打印出来
##usedmayi()
##print("[","蚂蚁类型号","，蚂蚁重量","，蚂蚁体积","，车型组合")
##print(Usedmayi)

###################################函数测试区域##################################################################
#################################################################################################################




#################################################################################################################

#################################################################################################################

#################################################################################################################

#################################################################################################################
    

