# RTK程序文档

### 结构体定义

##### 观测值结构体：

`obs_t`该结构体包含：

- 观测值是否读取进来了，即是否有效 `valid`

- 观测历元时间  包含（MJD和CT两种）
- 观测卫星数  `n`
- 存储每颗卫星观测值的数组 `obsd[MAXSATNUM]`

内置了 `Print()`函数，可往控制台打印出每个历元的观测值

**关于怎么索引某颗卫星观测值的问题**：根据PRN来索引，0-31 GPS  32-93 BDS，使用之前用`valid`来筛选是否有效



`obsd_t ` 该结构体包含：

- 系统标识符  `char sys`
- 卫星PRN号  `int prn`
- 双频伪距、载波、多普勒、信噪比

内置了 `Print()` 函数，可供 `obs_t`中的 `Print()`调用



`obs_header`该结构体包含表示GPS/BDS 需读取观测值类型的字符串，用于后续指定读取的顺序

- `string obsType[MAXOBSSYSTEM]`

内置了计算观测值顺序的函数 `void calobsTypeSequence(int obsTypeSequence[MAXOBSSYSTEM][8])`  这里计算出来的顺序返回至传进来的数组。8是因为双频、而每个频点4个观测值。



##### 星历结构体：

`eph_t ` 该结构体存放某颗卫星某个历元的星历（目前GPS和BDS存放的是一样的，后续需要加以改进）

- 九参数+六个摄动项
- 钟偏、钟速、钟漂
- `toc`(就是星历上第一行的时间戳)
- 系统标识符和PRN
- 传播时间 `double ttr` 
- 周数 `int week`  和周内秒 `double sow`



`nav_header` 导航电文头文件结构体，目前为空



`nav_t` 目前为空



##### 计算结果结构体：

`sat_t` 存储每颗卫星的位置、速度的结构体

- 有效位 `valid`

- 系统号`sys` 和 `PRN`
- 卫星信号发射时刻  `MJD transmitTime`
- 传播时间 ttr
- 三维速度、位置、姿态
- 钟差

存在的问题：不太确定ttr是否要保存在这个结构体里面，感觉对于ttr来说，应该保存在观测值里。

用处：

1. 用于保存计算卫星轨道位置、卫星测速的结果





### 读取Rinex文件

##### 观测值读取

**O文件头读取函数**  `void readObsHeader(std::fstream& fs, obs_header &obsHeader);` 

Param：1.文件流对象；2.用于保存头文件读取内容的`obs_header`类型的变量

Function：我们要能自动正确识别不同类型观测值在文件中的顺序，必须要读取#SYS OBS TYPE那一行头文件。读取之后保存了这些字符串，并且在函数内通过str的查找识别出CIC...的顺序，并将顺序返回至 `extern int OBSTYPESEQUENCE[MAXOBSSYSTEM][8];` 中，该变量中的顺序如"1 2 3 4 5 6 7 8"



O**文件内容读取函数** `void readObs(std::fstream& fs, obs_t &obs);`

Param：1.文件流对象；2.用于保存某一历元观测值数据的变量

Function: 逐历元读取观测值，并将其保存到`obs_t` 类型的变量中



2. 卫星轨道位置计算



3. 



### 构建函数模型

函数模型要输出的是矩阵B、$ L = o-c$  、待估参数阵x



### 最小二乘解算

设最多参与解算的卫星数为50颗

矩阵：

B  50*5

x  5*1    【x，y，z，dt_gps,  dt_bds】

l  o-c    伪距观测值-各项改正  -  两点计算出来的几何距离

P 单位阵

默认全0，有的才给他赋值。

x = $B^TPB^{-1}B^TPL$  



要有该历元的观测值obst(这里面含有改正)，卫星位置sat[]



### RTK函数模型构建

分别计算测站和流动站接收到的卫星位置，然后分别代入到B阵里面计算。



TO FIX 2021/5/30

NovaTel 0406 7190数据，输出大概8百多秒就出现问题了。

**是不是基准星发生了改变？，试着输出基准星**---应该不是

发生在sod=27804历元，l[3,14,25,36]达到10e6量级。

发现是8号星

**已解决**



浮点解大部分历元在0.1m之内，有少数几个历元最大到了4m左右。

加上高度角截止条件后，卫星数的确变少了，但是没看到浮点解精度变高

To fix 2021/5/31

浮点解没问题了，固定解还要改善，要迭代，把最后一次的改正数xyz和N放入lambda函数中。



##### 如何得到固定解：

首先我们已知流动站和基准站的概略坐标，认为基准站坐标是精确的，那么通过构造RTK模型，得到的最小二乘解里有关于流动站坐标的改正数，将其加到概略坐标上，再进行RTK模型的解算，直至某次xyz改正数小于给定的阈值退出，取最后一次的解算模糊度和方差作为lambda函数的输入，通过lambda方法固定整数模糊度，再通过式子代回去计算固定解的xyz改正数。