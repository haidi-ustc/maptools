# Materials Project Tool (maptool)
```
+--------------------------------------------------------------------+
|                     _   _   _   _   _   _   _                      |
|                    / \ / \ / \ / \ / \ / \ / \                     |
|                   ( m | a | p | t | o | o | l )                    |
|                    \_/ \_/ \_/ \_/ \_/ \_/ \_/                     |
|                                                                    |
|                     VERSION 0.1.1 (2018.10.29)                      |
|                                                                    |
|                   Written by Wang haidi(汪海迪)                    |
|                 URL  https://github.com/haidi-ustc                 |
|                 Bug report:(ustchaidi@gmail.com)                   |
+--------------------------------------------------------------------+
```


<html>


 <p><a href="https://github.com/haidi-ustc/image/blob/master/maptool_small.png?raw=true" target="_blank" rel="noopener noreferrer"><img src="https://github.com/haidi-ustc/image/raw/master/maptool_small.png?raw=true" height="300" style="max-width:100%;"></a>
<a href="https://github.com/haidi-ustc/image/blob/master/maptool_QR.jpg?raw=true" target="_blank" rel="noopener noreferrer"><img src="https://github.com/haidi-ustc/image/raw/master/maptool_QR.jpg?raw=true" height="300" style="max-width:100%;"></a></p>
</html>




> #### 使用maptool必须守以下license条款：
> * 用户可以自由下载，并免费使用 maptool 
> * 用户可以在此软件基础上针对自己的需求进行修改，但禁止传播/重新发布
> * 用户在作者未经授权情况下，禁止任何使用该软件进行的商业行为

# Installation
- install anaconda

联网安装 anaconda3（或者直接去anaconda官网下载linux x86_64版本）
    
```
wget https://repo.continuum.io/archive/Anaconda3-5.0.0-Linux-x86_64.sh
```
最后提示是否将环境变量加入```~/.bashrc```选择yes即可，当然别忘了执行
  
```
source ~/.bashrc
```
  
- 安装 atomate
```
pip install atomate
```
 
- 安装 maptool
```
python setup.py install 
```
- 设置 API_KEY



 用个人邮箱在materialsproject注册，点击右上角的```dashboad```即可获取API_KEY(内容类似下面的字符串)，然后用自己的API_KEY代替下面字符串```'fes20832fa9xfda28wa'```
 ```
 export MAPI_KEY='fes20832fa9xfda28wa'
 ```
 将上面内容写入```~/.bashrc```文件中，然后重启xshell或者执行
 ```
 source ~/.bashrc 
 ```
  
- 设置 VASP赝势 (可选步骤)

假定当前路径为/home/haidi/soft,在当前路径的vasp_pp文件中 将赝势文件按照下面方式解压存放（这里只显示部分文件）
 
```
 vasp_pp/
|-- potpaw_LDA
|   |-- Ac
|   |-- Ag
|   |-- Ag_GW
|   |-- Ag_pv
|   |-- data_base
|   `-- potpaw_LDA.54.tar.gz
`-- potpaw_PBE
    |-- Ac
    |-- Ag
    |-- Ag_GW
    |-- Ag_pv
    |-- Ag_sv_GW
    |-- Al
    `-- potpaw_PBE.5.4.tar.gz
```

```
cd /home/haidi/soft
pmg config -p /home/haidi/vasp_pp /home/haidi/.pymatgen
pmg config --add PMG_VASP_PSP_DIR /home/haidi/.pymatgen
```

完成上面操作之后执行
```
pmg potcar -s Si C
```
如果产生SiC的POTCAR，则表示赝势配置正常

- 测试

安装完成后执行
```
mpt -i
```
如果输出如下类似信息即可表示安装成功：
```
 True
Maptool
--------

Version: 0.1.0
Path:    /home/vasp/.local/lib/python3.6/site-packages/maptool-0.1.0-py3.6.egg/maptool
Date:    May 16, 2018

Python version=3.6.2 |Continuum Analytics, Inc.| (default, Jul 20 2017, 13:51:32) 
[GCC 4.4.7 20120313 (Red Hat 4.4.7-1)]

   pymongo      3.6.1   /home/vasp/.pyenv/versions/anaconda3-4.3.0/envs/mpt_a3/lib/python3.6/site-packages/pymongo
     numpy     1.14.3   /home/vasp/.pyenv/versions/anaconda3-4.3.0/envs/mpt_a3/lib/python3.6/site-packages/numpy
     scipy      1.0.1   /home/vasp/.pyenv/versions/anaconda3-4.3.0/envs/mpt_a3/lib/python3.6/site-packages/scipy
    mayavi            Not Found
matplotlib      2.2.2   /home/vasp/.pyenv/versions/anaconda3-4.3.0/envs/mpt_a3/lib/python3.6/site-packages/matplotlib
      tqdm     4.23.1   /home/vasp/.pyenv/versions/anaconda3-4.3.0/envs/mpt_a3/lib/python3.6/site-packages/tqdm
    future     0.16.0   /home/vasp/.pyenv/versions/anaconda3-4.3.0/envs/mpt_a3/lib/python3.6/site-packages/future
      nose      1.3.7   /home/vasp/.pyenv/versions/anaconda3-4.3.0/envs/mpt_a3/lib/python3.6/site-packages/nose
  coverage      4.5.1   /home/vasp/.pyenv/versions/anaconda3-4.3.0/envs/mpt_a3/lib/python3.6/site-packages/coverage
    spglib     1.10.3   /home/vasp/.pyenv/versions/anaconda3-4.3.0/envs/mpt_a3/lib/python3.6/site-packages/spglib
    pyhull   2015.2.0   /home/vasp/.pyenv/versions/anaconda3-4.3.0/envs/mpt_a3/lib/python3.6/site-packages/pyhull
  pymatgen  2018.9.12   /home/vasp/.pyenv/versions/anaconda3-4.3.0/envs/mpt_a3/lib/python3.6/site-packages/pymatgen
      qmpy            Not Found
       ase     3.16.0   /home/vasp/.pyenv/versions/anaconda3-4.3.0/envs/mpt_a3/lib/python3.6/site-packages/ase
 ```

> 注： matpool的可执行文件名不是```maptool```而是 ```mpt```！！！，所以运行的时候不带参数时，命令行只需输入：
```
mpt
```
![maptool](https://github.com/haidi-ustc/image/blob/master/maptool-0.1.1.png?raw=true)

# Maptool V.S. Vaspkit
下面表格对maptool和vaspkit进行了简单对比，二者各有优缺点，maptool最大的优点在于在线数据查询以及快速的数据可视化，但是数据查询需要联网（对于国内有些超算可能不大方便），此外maptool安装过程需要联网以简化安装步骤（不联网理论上可以安装，但是相对比较复杂，涉及到非常多的python函数库），相比之下vaspkit由于采用Fortran编写，所以安装比较简单，但也导致其开发难度较大，很难实现高级功能。

PK | maptool | vaspkit
---|---| ---
开发语言 | python | Fortran
支持环境 | linux  | linux
数据库查询| 支持   | 不支持
数据可视化 | 支持  | 不支持
支持软件 | VASP | VASP
能带高对称点 | 任意结构  | 特定结构
高通量输出  | 预留接口   | 无
运行速度  | 稍慢 (无关紧要)   | 快


# Menu

- [======================== structural operation ========================](#structural-operation) 
  - [x] [a1 >>> random operation](#random-operation)
  - [x] [a2 >>> covert operation](#convert-operation)
  - [x] [a3 >>> build operation  ](#build-operation)
  - [x] [a4 >>> cleave operation ](#cleave-operation)
  - [x] [a5 >>> strain operation ](#strain-operation)
  - [x] [a6 >>> TwoD operation](#TwoD-operation)
- [======================== structural analysis =========================](#structural-analysis)
  - [x] [b1 >>> structure symmetry ](#structure-symmetry)
  - [x] [b2 >>> structure finger print](#structure-finger-print)
  - [x] [b3 >>> structure difference](#structure-difference)
  - [x] [b4 >>> get primitive cell](#get-primitive-cell)
  - [x] [b5 >>> get conventional cell](#get-conventional-cell)
  - [x] [b6 >>> get XRD pattern](#get-XRD-pattern)
- [========================= vasp in/out tools ==========================](#vasp-in/out-tools)
  - [x] [c1 >>> prepare input files](#prepare-input-files)
  - [x] [c2 >>> analysis output files](#analysis-output-files)
- [======================= vasp auto calculation ========================](#vasp-auto-calculation)
  - [ ] [d1 >>> optimize structure](#optimize-structure)
  - [ ] [d2 >>> calculate band structure](#calculate-band-structure)
  - [ ] [d3 >>> calculate band structure HSE06](#calculate-band-structure-HSE06)
  - [ ] [d4 >>> calculate dos](#calculate-dos)
  - [ ] [d5 >>> calculate dos by HSE06](#calculate-dos-by-HSE06)
  - [ ] [d6 >>> calculate elastic properties](#calculate-elastic-properties)
  - [ ] [d7 >>> calculate phonon](#calculate-phonon)
  - [ ] [d8 >>> execute MD simulation](#execute-MD-simulation)
- [=========================  online retrieve =========================](#online-retrieve)
  - [x] [e1 >>> get band/dos by mp-ID](#get-band/dos-by-mp-ID)
  - [x] [e2 >>> get structure from materialsproject database](#get-structure-from-materialsprojecet-database)
  - [x] [e3 >>> get properties by mp-ID](#get-properties-by-mp-ID)


> -  maptool目前版本只支持菜单选择模式（后期可能会开发命令行参数和图形界面），所以在使用过程中只需要对照相应的菜单，输入前面所指示的数字，字母或者数字字母的组合。对于有些需要非常多参数的操作，要求使用者编写相应的输入文件，以见命令行参数的输入个数。
> -  在读取以```POSCAR```格式的晶体结构文件的过程中，务必保证文件名以POSCAR开头或者以POSCAR结尾，例如：POSCAR-2和a_POSCAR符合要求，而my-POSCAR-1不合理。
> -  所有--------------->>下面一行的字符均需要用户输入。
> - 以上打钩选项表示该功能可用。
> - 用户屏幕输出和手册显示内容若有不一致可能是因为输入文件不用导致。
----


# structure #

maptool 内部所有的结构以[pymatgen](www.materialsproject.org)中的```Structure```类为基础进行操作。读取过程分成周期性和非周期性结构，需要手动指定结构是否具有周期。

周期性结构：

```
your choice ?
1 >>> crystal
2 >>> molecule/cluster
--------------->>
1
your choice ?
0 >>> specific the file name!
1 >>> POSCAR
2 >>> CONTCAR
3 >>> CHGCAR
4 >>> LOCPOT
5 >>> maptool.cif
6 >>> maptool.xsf
7 >>> maptool.nc
8 >>> maptool.yaml
9 >>> maptool.json
--------------->>
1
```
上述操作表明结构将从```POSCAR```的文件中读取


非周期性结构：

```
your choice ?
1 >>> crystal
2 >>> molecule/cluster
--------------->>
2
your choice ?
0 >>> specific the file name!
1 >>> maptool.xyz
2 >>> maptool.mol
3 >>> maptool.nc
4 >>> maptool.yaml
5 >>> maptool.json
--------------->>
0
input the full file name
--------------->>
a.xyz
```
上述操作表明结构将从```a.xyz```的文件中读取

如果上述第二个选择的时候输入1，则表示结构将从```maptool.xyz```文件中读取，剩下的以此类推。

----

# structural operation

结构操作模块，主要功能包括对晶体和分子结构的随机产操作，转换操作，构建操作，应变操作和二维结构操作。

## random operation

- ### random structure generating

 产生随机结构
 
 ##### 命令行参数：
 ```
your choice ?
--------------->>
a1
your choice ?
1 >>> random structure generating
2 >>> random perturbation for atom index
3 >>> random disturbing for lattice matrix
4 >>> random disturbing for atom position
--------------->>
1
input the formula of structure like this: 
(Fe3O4)4
or input the formula and spacegroup like this 
Si4O8  20
--------------->>
(Fe3O4)4
```
表示产生包括4个单元的```$Fe_3O_4$```的结构，即```$Fe_{12}O_{16}$```的随机结构，或者
```
your choice ?
--------------->>
a1
your choice ?
1 >>> random structure generating
2 >>> random perturbation for atom index
3 >>> random disturbing for lattice matrix
4 >>> random disturbing for atom position
--------------->>
1
input the formula of structure like this: 
(Fe3O4)4
or input the formula and spacegroup like this 
Si4O8  20
--------------->>
Si4O8  20
```
表示用第20号空间群来产生```$ Si_{4}O_8 $```结构。

##### 输出文件：

random_spg_群号.vasp

- ### random perturbation for atom index


随机置换两个原子

##### 命令行参数：

```
your choice ?
--------------->>
a1
your choice ?
1 >>> random structure generating
2 >>> random perturbation for atom index
3 >>> random disturbing for lattice matrix
4 >>> random disturbing for atom position
--------------->>
2
your choice ?
1 >>> crystal
2 >>> molecule/cluster
--------------->>
1
your choice ?
0 >>> specific the file name!
1 >>> POSCAR
2 >>> CONTCAR
3 >>> CHGCAR
4 >>> LOCPOT
5 >>> maptool.cif
6 >>> maptool.xsf
7 >>> maptool.nc
8 >>> maptool.yaml
9 >>> maptool.json
--------------->>
1
```

##### 输出文件：

swap_原子1序号_原子2序号.vasp


- ### random disturbing for lattice matrix
随机扰动晶胞参数```$\alpha,\beta,\gamma,a,b,c$```,表现为晶胞矩阵的变化


##### 命令行参数：
```
your choice ?
--------------->>
a1
your choice ?
1 >>> random structure generating
2 >>> random perturbation for atom index
3 >>> random disturbing for lattice matrix
4 >>> random disturbing for atom position
--------------->>
3
your choice ?
0 >>> specific the file name!
1 >>> POSCAR
2 >>> CONTCAR
3 >>> CHGCAR
4 >>> LOCPOT
5 >>> maptool.cif
6 >>> maptool.xsf
7 >>> maptool.nc
8 >>> maptool.yaml
9 >>> maptool.json
--------------->>
1
input the maximum displacement and with fix diagonal 
or non-diagonal element like this:  0.01 T F
it means maximum displacement is 0.01 
the diagonal element will be fixed,
while random disturbing will be add to non-digaonal element
--------------->>
0.01 T F
```
表示最大‘位移’为0.01，第一个T表示固定对角项，第二个F表示变化非对角项。

##### 输出文件：

random.vasp
 
- ### random disturbing for atom position

随机移动一个原子的位置



##### 命令行参数：
```
your choice ?
--------------->>
a1
your choice ?
1 >>> random structure generating
2 >>> random perturbation for atom index
3 >>> random disturbing for lattice matrix
4 >>> random disturbing for atom position
--------------->>
4
input the maximum displacement(<0.25 in Angstrom)
--------------->>
0.025
your choice ?
1 >>> crystal
2 >>> molecule/cluster
--------------->>
1
your choice ?
0 >>> specific the file name!
1 >>> POSCAR
2 >>> CONTCAR
3 >>> CHGCAR
4 >>> LOCPOT
5 >>> maptool.cif
6 >>> maptool.xsf
7 >>> maptool.nc
8 >>> maptool.yaml
9 >>> maptool.json
--------------->>
1
```
表示原子移动的真实位移量为0.025 Ang

##### 输出文件：

random.vasp
 
 
## convert operation

结构转换模块，一般用于POSCAR为其他格式，或者其他格式转换为POSAR
> 注意： 输出文件的文件名必须包括文件格式后缀

## build operation

- ### build supercell
扩建超胞结构


##### 命令行参数：
```
your choice ?
--------------->>
a3
your choice ?
1 >>> build supercell
2 >>> build nanotube
3 >>> build absorption configuration
1
your choice ?
0 >>> specific the file name!
1 >>> POSCAR
2 >>> CONTCAR
3 >>> CHGCAR
4 >>> LOCPOT
5 >>> maptool.cif
6 >>> maptool.xsf
7 >>> maptool.nc
8 >>> maptool.yaml
9 >>> maptool.json
--------------->>
1
--------------->>

Several options are available:

a. A full 3x3 scaling matrix defining the linear combination
   the old lattice vectors. E.g., 2 1 0  0 1 0  0 0 3
   generates a new structure with lattice vectors a' =
   2a + b, b' = 3b, c' = c where a, b, and c are the lattice
   vectors of the original structure.
b. An sequence of three scaling factors. E.g., 2 1 1
   specifies that the supercell should have dimensions 2a x b x
   c.
c. A number, which simply scales all lattice vectors by the
   same factor.
        
--------------->>
2 1 0  0 1 0  0 0 3
```
表示```$ a' = 2a + b, b' = 3b, c'= c$``` 这里的```a,b,c```表示原始晶胞的三个矢量，```a',b',c'```表示新结构的晶胞矢量，或者
##### 命令行参数：
```
your choice ?
--------------->>
a3
your choice ?
1 >>> build supercell
2 >>> build nanotube
3 >>> build absorption configuration
1
your choice ?
0 >>> specific the file name!
1 >>> POSCAR
2 >>> CONTCAR
3 >>> CHGCAR
4 >>> LOCPOT
5 >>> maptool.cif
6 >>> maptool.xsf
7 >>> maptool.nc
8 >>> maptool.yaml
9 >>> maptool.json
--------------->>
1
--------------->>

Several options are available:

a. A full 3x3 scaling matrix defining the linear combination
   the old lattice vectors. E.g., 2 1 0  0 1 0  0 0 3
   generates a new structure with lattice vectors a' =
   2a + b, b' = 3b, c' = c where a, b, and c are the lattice
   vectors of the original structure.
b. An sequence of three scaling factors. E.g., 2 1 1
   specifies that the supercell should have dimensions 2a x b x
   c.
c. A number, which simply scales all lattice vectors by the
   same factor.
        
--------------->>
 2 1 1
```
表示产生```2x2x1```的超胞，或者

##### 命令行参数：
```
your choice ?
--------------->>
a3
your choice ?
1 >>> build supercell
2 >>> build nanotube
3 >>> build absorption configuration
1
your choice ?
0 >>> specific the file name!
1 >>> POSCAR
2 >>> CONTCAR
3 >>> CHGCAR
4 >>> LOCPOT
5 >>> maptool.cif
6 >>> maptool.xsf
7 >>> maptool.nc
8 >>> maptool.yaml
9 >>> maptool.json
--------------->>
1
--------------->>

Several options are available:

a. A full 3x3 scaling matrix defining the linear combination
   the old lattice vectors. E.g., 2 1 0  0 1 0  0 0 3
   generates a new structure with lattice vectors a' =
   2a + b, b' = 3b, c' = c where a, b, and c are the lattice
   vectors of the original structure.
b. An sequence of three scaling factors. E.g., 2 1 1
   specifies that the supercell should have dimensions 2a x b x
   c.
c. A number, which simply scales all lattice vectors by the
   same factor.
        
--------------->>
3
```
表示产生```3x3x3```的超胞

##### 输出文件：

supercell.vasp

- ###  build nanotube

暂时没有实现

- ### build absorption configuration

自动构建分子在衬底上的吸附模型，程序自动计算吸附位点，给出各种可能的吸附构型，以POSCAR格式出文件

命令行参数格式：
```
your choice ?
--------------->>
a3
your choice ?
1 >>> build supercell
2 >>> build nanotube
3 >>> build absorption configuration
3
your choice ?
1 >>> read slab from file
2 >>> build slab by bulk
--------------->>
2
your choice ?
0 >>> specific the file name!
1 >>> POSCAR
2 >>> CONTCAR
3 >>> CHGCAR
4 >>> LOCPOT
5 >>> maptool.cif
6 >>> maptool.xsf
7 >>> maptool.nc
8 >>> maptool.yaml
9 >>> maptool.json
--------------->>
1
your choice ?
0 >>> specific the file name!
1 >>> maptool.xyz
2 >>> maptool.mol
3 >>> maptool.nc
4 >>> maptool.yaml
5 >>> maptool.json
--------------->>
0
input the full file name
--------------->>
a.xyz
```
表示衬底从```POSCAR```文件读取，吸附分子文件从```a.xyz```文件读取

##### 输出文件：

adsorb2_0_1-15.vasp （表示吸附在结构201面上的第15个结构 ）

或者

```
your choice ?
--------------->>
a3
your choice ?
1 >>> build supercell
2 >>> build nanotube
3 >>> build absorption configuration
--------------->>
3
your choice ?
1 >>> read slab from file
2 >>> build slab by bulk
--------------->>
1
your choice ?
0 >>> specific the file name!
1 >>> POSCAR
2 >>> CONTCAR
3 >>> CHGCAR
4 >>> LOCPOT
5 >>> maptool.cif
6 >>> maptool.xsf
7 >>> maptool.nc
8 >>> maptool.yaml
9 >>> maptool.json
--------------->>
0
input the full file name
--------------->>
POSCAR_slab
your choice ?
0 >>> specific the file name!
1 >>> maptool.xyz
2 >>> maptool.mol
3 >>> maptool.nc
4 >>> maptool.yaml
5 >>> maptool.json
--------------->>
0
input the full file name
--------------->>
a.xyz
```
表示衬底从POSCAR_slab文件中读取，吸附分子从a.xyz文件中读取。

##### 输出文件：

 adsorb-序号.vasp

> 上面方法产生的吸附构型并没有明确指定吸附的面，衬底是否扩包等条件，为了更加明确的控制产生的吸附构型，则可以用配置文件来指定。


以水分子吸附在Ni表面为例，来说明配置文件的编写规则。

##### 输入文件:

adsorb.cfg  h2o.xyz  Ni.POSCAR

其中配置文件的内容为：
```
method    = 2              # method 
crystal   = Ni.POSCAR      # bulk or slab file
molecule  = h2o.xyz        # molecule file
max_index = 2              # max miller index
min_vacum = 20             # minimum vacuum slab thickness
min_slab  = 8              # mininum slab thickness
repeat    = 3 3 1          # repeat unit for substrate
```
method  = 2 表示衬底从体相结构构建

method  = 1 表示衬底从slab文件读取

上述配置文件表示通过Ni.POSCAR文件读取晶体文件并且从概文件构建不同的晶面，从h2o.xyz文件读取水分子构型，最多晶面指数序列为2，最小真空层20 埃，衬底由``` $ 3X3X1 $ ```的超胞构成。

命令行运行参数：
```
your choice ?
--------------->>
a3
your choice ?
1 >>> build supercell
2 >>> build nanotube
3 >>> build absorption configuration
--------------->>
3
```
这里输入3之后程序会自动从adsorb.cfg文件去读取参数，构建吸附模型的POSCAR

## cleave operation

- ### cleave surface
晶体切面操作

##### 输入文件：
POSCAR

##### 命令行参数：
```
your choice ?
--------------->>
a4
your choice ?
0 >>> specific the file name!
1 >>> POSCAR
2 >>> CONTCAR
3 >>> CHGCAR
4 >>> LOCPOT
5 >>> maptool.cif
6 >>> maptool.xsf
7 >>> maptool.nc
8 >>> maptool.yaml
9 >>> maptool.json
--------------->>
1
your choice ?
1 >>> cleave surface
2 >>> cleave ball
3 >>> cleave ball shell
--------------->>
1
 input the miller index, minimum size in angstroms of layers containing atomssupercell
 and Minimize size in angstroms of layers containing vacuum like this:
 1 0 0 | 5 | 5
 it means miller index is [1,0,0]
 min_slab_size is 5 Ang 
 min_vacum_size is 5 Ang 
 or like this : 
 2 | 5 | 5
 it will generate all slab with miller index less than 2
--------------->>
1 0 0 | 5 | 5
```
这里的```1 0 0 | 5 | 5```表示产生100晶面，slab大小为5埃，真空层大小为5埃，注意用空格和竖线(|)分割数据
##### 输出文件：
1_0_0.vasp

也可以按下面格式输入

##### 输入文件：
POSCAR

##### 命令行参数：
```
your choice ?
--------------->>
a4
your choice ?
0 >>> specific the file name!
1 >>> POSCAR
2 >>> CONTCAR
3 >>> CHGCAR
4 >>> LOCPOT
5 >>> maptool.cif
6 >>> maptool.xsf
7 >>> maptool.nc
8 >>> maptool.yaml
9 >>> maptool.json
--------------->>
1
your choice ?
1 >>> cleave surface
2 >>> cleave ball
3 >>> cleave ball shell
--------------->>
1
 input the miller index, minimum size in angstroms of layers containing atomssupercell
 and Minimize size in angstroms of layers containing vacuum like this:
 1 0 0 | 5 | 5
 it means miller index is [1,0,0]
 min_slab_size is 5 Ang 
 min_vacum_size is 5 Ang 
 or like this : 
 2 | 5 | 5
 it will generate all slab with miller index less than 2
--------------->>
 2 | 5 | 5
```
这里的```2 | 5 | 5```表示产生所有晶面指数序列小于2的晶面，slab大小为5埃，真空层大小为5埃，注意用空格和竖线(|)分割数据

##### 输出文件：
h_k_l.vasp

- ### cleave ball

以某个原子为中心，截取指定半径的球体，转换为POSCAR格式


 ##### 输入文件：
 POSCAR

 ##### 命令行参数：
 
```
your choice ?
--------------->>
a4
your choice ?
0 >>> specific the file name!
1 >>> POSCAR
2 >>> CONTCAR
3 >>> CHGCAR
4 >>> LOCPOT
5 >>> maptool.cif
6 >>> maptool.xsf
7 >>> maptool.nc
8 >>> maptool.yaml
9 >>> maptool.json
--------------->>
1
your choice ?
1 >>> cleave surface
2 >>> cleave sphere cluster
3 >>> cleave shell structure
--------------->>
2
 input the center atom index, sphere radius and vacuum layer thickness
 1 3.5 15
 it means the sphere will be selected according to the 1st atom
 with the radius equals 5Ang, and vacuum layer thickness is 15 Ang
--------------->>
1 3.5 15
```
以1号原子为中心，截取球体，半径3.5A， 真空层15A

 ##### 输出文件:
 
 sphere.vasp 
 

- ### cleave ball shell

截取指定厚度的球壳模型，并转换为POSCAR格式

![NFE](https://github.com/haidi-ustc/image/blob/master/shell.png?raw=true)

 ##### 输入文件：
 POSCAR

 ##### 命令行参数：
 

```
your choice ?
--------------->>
a4
your choice ?
0 >>> specific the file name!
1 >>> POSCAR
2 >>> CONTCAR
3 >>> CHGCAR
4 >>> LOCPOT
5 >>> maptool.cif
6 >>> maptool.xsf
7 >>> maptool.nc
8 >>> maptool.yaml
9 >>> maptool.json
--------------->>
1
your choice ?
1 >>> cleave surface
2 >>> cleave sphere cluster
3 >>> cleave shell structure
--------------->>
3
 input the center atom index, start radius, shell thickness and
 vacuum layer thickness
 1 5 10  15
 it means the ball shell will be selected according to the 1st atom
 with the 5< r <15Ang, and vacuum layer thickness is 15 Ang
--------------->>
1 10 1.2 15
```
以1号原子为中心，截取球壳，空壳半径10A， 球壳厚度1.2A， 真空层15A


 ##### 输出文件:
 
shell.vasp 
 

## strain operation

对晶体结构施加应变操作，得到一系列应变后的结构

##### 输入文件：
POSCAR

##### 命令行参数：
```
your choice ?
--------------->>
a5
your choice ?
0 >>> specific the file name!
1 >>> POSCAR
2 >>> CONTCAR
3 >>> CHGCAR
4 >>> LOCPOT
5 >>> maptool.cif
6 >>> maptool.xsf
7 >>> maptool.nc
8 >>> maptool.yaml
9 >>> maptool.json
--------------->>
1
input the strain component like :
0.01
it means aplly a strain of 1% along all directions
or
0.01 0.0 0.0
it means apply a strain of 1% along the x direction
or
0.01:0.03:5 0.0 0.0
it means to devide strain range into 5 parts
--------------->>
0.01
```
表示对原有结构每个方向施加1%的应变(晶体格子整体放大1%)

##### 输出文件：
index_strain.dat    

strain_0.vasp

或者按照下面格式输入
##### 输入文件：
POSCAR

##### 命令行参数：
```
your choice ?
--------------->>
a5
your choice ?
0 >>> specific the file name!
1 >>> POSCAR
2 >>> CONTCAR
3 >>> CHGCAR
4 >>> LOCPOT
5 >>> maptool.cif
6 >>> maptool.xsf
7 >>> maptool.nc
8 >>> maptool.yaml
9 >>> maptool.json
--------------->>
1
input the strain component like :
0.01
it means aplly a strain of 1% along all directions
or
0.01 0.0 0.0
it means apply a strain of 1% along the x direction
or
0.01:0.03:5 0.0 0.0
it means to devide strain range into 5 parts
--------------->>
0.01 0.0 0.0
```
表示对原有结构```$x$```方向施加1%的应变

##### 输出文件：
index_strain.dat   

strain_0.vasp


或者产生一系列应变的结构

##### 输入文件：
POSCAR

##### 命令行参数：

```
your choice ?
--------------->>
a5
your choice ?
0 >>> specific the file name!
1 >>> POSCAR
2 >>> CONTCAR
3 >>> CHGCAR
4 >>> LOCPOT
5 >>> maptool.cif
6 >>> maptool.xsf
7 >>> maptool.nc
8 >>> maptool.yaml
9 >>> maptool.json
--------------->>
1
input the strain component like :
0.01
it means aplly a strain of 1% along all directions
or
0.01 0.0 0.0
it means apply a strain of 1% along the x direction
or
0.01:0.03:5 0.0 0.0
it means to devide strain range into 5 parts
--------------->>
0.01:0.03:5 0.0 0.0 
```


表示对原有结构```$x$```方向施加1%到3%的应变，分割为5个点

##### 输出文件：
index_strain.dat   

strain_0.vasp

strain_1.vasp

strain_2.vasp  

strain_3.vasp  

strain_4.vasp 

每个结构对于的应变```$(r-r_0)/r_0$```写到index_strain.dat文件中，其
内容格式大致如下：

```
  0  0.0100  0.0000  0.0000
  1  0.0150  0.0000  0.0000
  2  0.0200  0.0000  0.0000
  3  0.0250  0.0000  0.0000
  4  0.0300  0.0000  0.0000
 ```


## TwoD operation
2D材料相关操作模块
- ### build rippled structure
产生褶皱2D结构,参考文献 [Anisotropic Ripple Deformation in Phosphorene](https://pubs.acs.org/doi/abs/10.1021/acs.jpclett.5b00522)

![ripple phosphorene](https://pubs.acs.org/appl/literatum/publisher/achs/journals/content/jpclcd/2015/jpclcd.2015.6.issue-9/acs.jpclett.5b00522/20150501/images/medium/jz-2015-00522b_0006.gif)


##### 输入文件：
POSCAR，要求单胞为正交结构

##### 命令行参数：
```
your choice ?
--------------->>
a6
your choice ?
1 >>> build rippled structure
2 >>> build multi-layered structure
3 >>> split multi-layered structure
4 >>> resize vacuum layer
5 >>> center atomic-layer along z direction
6 >>> apply strain along different direction
7 >>> constrain atom in specific range
8 >>> get a substrate for 2D material (online!!!)
--------------->>
1
your choice ?
0 >>> specific the file name!
1 >>> POSCAR
2 >>> CONTCAR
3 >>> CHGCAR
4 >>> LOCPOT
5 >>> maptool.cif
6 >>> maptool.xsf
7 >>> maptool.nc
8 >>> maptool.yaml
9 >>> maptool.json
--------------->>
1
input the supercell scaling factor
for x direction can be: 10 1 1
for y direction can be: 1 10 1
--------------->>
10 1 1
input the strain range
example: 0.02:0.1:10 
--------------->>
0.02:0.1:10
input the index of atom need to be fixed
example: 1 10 11 20 
0 means fix the atom automatically
--------------->>
0
```
为了构建ripple结构，首先需要沿着某一方向建立超胞，这里```10 1 1```表示建立```10x1x1```的超胞，```0.02:0.1:10```表示沿着```y```方向施加2%到10%的应变，分割为10份。命令行输入的最后告诉程序边缘的那些原子需要固定，0表示由程序自动确定边缘需要固定的原子。

>  这里之说以要施加应变，因为对材料施加应变之后才有可能使其出现褶皱，具体参考上面文献。

##### 输出文件：

 应变分数_wo.vasp   #只施加了应变的结构
 
 应变分数_w.vasp    #施加了应变，发生褶皱的结构
 
 

- ### build multi-layered structure
由单层结构构建多层结构


##### 输入文件：
POSCAR，要求c轴垂直于x-y平面

##### 命令行参数：

```
your choice ?
--------------->>
a6
your choice ?
1 >>> build rippled structure
2 >>> build multi-layered structure
3 >>> split multi-layered structure
4 >>> resize vacuum layer
5 >>> center atomic-layer along z direction
6 >>> apply strain along different direction
7 >>> constrain atom in specific range
8 >>> get a substrate for 2D material (online!!!)
--------------->>
2
your choice ?
0 >>> specific the file name!
1 >>> POSCAR
2 >>> CONTCAR
3 >>> CHGCAR
4 >>> LOCPOT
5 >>> maptool.cif
6 >>> maptool.xsf
7 >>> maptool.nc
8 >>> maptool.yaml
9 >>> maptool.json
--------------->>
1
input the number of layers
--------------->>
2
input the layer distance
--------------->>
3.0
```
这里从POSCAR读取结构，然后构建其双层结构，层间距为3.0埃。

##### 输出文件：
layer_2.vasp


- ### split multi-layered structure
剥离多层结构中的一层或者多层，用于计算剥离能,参考文献 [Obtaining Two-Dimensional Electron Gas in Free Space without Resorting to Electron Doping: An Electride Based Design](https://pubs.acs.org/doi/pdf/10.1021/ja5065125)

![NFE](https://github.com/haidi-ustc/image/blob/master/zhaost_cleave.png?raw=true)


##### 输入文件：
POSCAR，要求z方向为真空层方向，c轴至少50埃

##### 命令行参数：
```
your choice ?
--------------->>
a6
========================2D structure operation========================
your choice ?
1 >>> build rippled structure
2 >>> build multi-layered structure
3 >>> split multi-layered structure
4 >>> resize vacuum layer
5 >>> center atomic-layer along z direction
6 >>> apply strain along different direction
7 >>> constrain atom in specific range
8 >>> get a substrate for 2D material (online!!!)
--------------->>
3
your choice ?
0 >>> specific the file name!
1 >>> POSCAR
2 >>> CONTCAR
3 >>> CHGCAR
4 >>> LOCPOT
5 >>> maptool.cif
6 >>> maptool.xsf
7 >>> maptool.nc
8 >>> maptool.yaml
9 >>> maptool.json
--------------->>
1

    input data according to tips

     select atoms by following ways:
     1. atomic index in POSCAR
        i.e. : 1 2 4-8 10 12-30
        i.e. : 1 2 4 8 10 
     2. atomic label
        i.e. : Si O
     3. atomic position
        i.e. : 0 0.5 | 0.2 0.4 | 0.3 0.7
        this means atoms with 0<x<0.5, 
        0.2<y<0.4 and 0.3<z<0.7 will be seleted
        or just specific the z coordinates,
        i.e. : ||0.3 0.7
        
--------------->>
1 2 4-8 10 12-30          
input the splitting distance, 10 Ang is enough!
--------------->>
10
numbers of splitting site, 50 sites are enough!
--------------->>
50
```
上述操作用于读取POSCAR文件，选择其中1,2,4到8,10,12到30号原子作为移动目标，最终移动10埃，分割为50步进行移动。

> #### maptool存在三种方式选择特定的原子：

>   - 原子序号选择，可以用 -表示连续的范围，用空格分开
>   - 元素符号，用空格分开
>   - 原子位置，用 | 分割，分别表示x，y和z的范围，如果为空，则表示该方向不受限制

##### 输出文件：
序号.vasp

split.dat

其中，split.dat文件格式如下：
```
#       index  distance/Ang
           0     0.000000
           1     1.750000
           2     3.500000
           3     2.750000
           4     2.000000
```
第一列表示序号，第二列表示剥离距离


- ### resize vacuum layer
改变真空层的厚度
##### 输入文件：
POSCAR，要求z方向为真空层方向
##### 命令行参数：
```
your choice ?
--------------->>
a6
your choice ?
1 >>> build rippled structure
2 >>> build multi-layered structure
3 >>> split multi-layered structure
4 >>> resize vacuum layer
5 >>> center atomic-layer along z direction
6 >>> apply strain along different direction
7 >>> constrain atom in specific range
8 >>> get a substrate for 2D material (online!!!)
--------------->>
4
your choice ?
0 >>> specific the file name!
1 >>> POSCAR
2 >>> CONTCAR
3 >>> CHGCAR
4 >>> LOCPOT
5 >>> maptool.cif
6 >>> maptool.xsf
7 >>> maptool.nc
8 >>> maptool.yaml
9 >>> maptool.json
--------------->>
1
current vacuum layer thickness is  4.759 Ang
input the new value of vacuum layer thickness
--------------->>
20
```
读取POSCAR文件，讲原有的真空层厚度4.759埃改成20埃。
##### 输出文件：
new_vacuum.vasp



- ### center atomic-layer along z direction
移动结构到到中心位置0.5 
##### 输入文件：
POSCAR，要求z方向为真空层方向
##### 命令行参数：
```
your choice ?
--------------->>
a6
your choice ?
1 >>> build rippled structure
2 >>> build multi-layered structure
3 >>> split multi-layered structure
4 >>> resize vacuum layer
5 >>> center atomic-layer along z direction
6 >>> apply strain along different direction
7 >>> constrain atom in specific range
8 >>> get a substrate for 2D material (online!!!)
--------------->>
5
your choice ?
0 >>> specific the file name!
1 >>> POSCAR
2 >>> CONTCAR
3 >>> CHGCAR
4 >>> LOCPOT
5 >>> maptool.cif
6 >>> maptool.xsf
7 >>> maptool.nc
8 >>> maptool.yaml
9 >>> maptool.json
--------------->>
1
```
读取POSCAR文件，将其移动到z-方向的中心位置。
##### 输出文件：
z-center.vasp


- ### apply strain along different direction
对2D结构在面内施加不同方向的应力，参考文献 [δ-Phosphorene: a two dimensional material with a highly negative Poisson's ratio](http://pubs.rsc.org/en/content/articlelanding/2017/nr/c6nr08550d/unauth#!divAbstract)

![delta-P](https://github.com/haidi-ustc/image/blob/master/hd_strain.png?raw=true)

 
##### 输入文件：
POSCAR，要求z方向为真空层方向，正交结构

##### 命令行参数：
```
your choice ?
--------------->>
a6
your choice ?
1 >>> build rippled structure
2 >>> build multi-layered structure
3 >>> split multi-layered structure
4 >>> resize vacuum layer
5 >>> center atomic-layer along z direction
6 >>> apply strain along different direction
7 >>> constrain atom in specific range
8 >>> get a substrate for 2D material (online!!!)
--------------->>
6
your choice ?
0 >>> specific the file name!
1 >>> POSCAR
2 >>> CONTCAR
3 >>> CHGCAR
4 >>> LOCPOT
5 >>> maptool.cif
6 >>> maptool.xsf
7 >>> maptool.nc
8 >>> maptool.yaml
9 >>> maptool.json
--------------->>
1
input the elastic of material by order : C11 C12 C22 C66
--------------->>
88.64 -23.71 149.21 24.50      
input applied force: e.x. 1.0 GPa nm
--------------->>
1.0
+-----
```
从POSCAR读取结构，要求输入弹性常数，顺序为```$C_{11}, C_{12} ,C_{22},C_{66}$```，最后需要输入每个方向施加的应力，具体定义参数上述文献。

##### 输出文件：
0.vasp - 36.vasp

即沿着360度的方向施加压力，产生36个文件

- ### constrain atom in specific range
固定原子坐标，用于限制性优化

##### 输入文件：
POSCAR

##### 命令行参数:
```
your choice ?
--------------->>
a6
your choice ?
1 >>> build rippled structure
2 >>> build multi-layered structure
3 >>> split multi-layered structure
4 >>> resize vacuum layer
5 >>> center atomic-layer along z direction
6 >>> apply strain along different direction
7 >>> constrain atom in specific range
8 >>> get a substrate for 2D material (online!!!)
--------------->>
7
your choice ?
0 >>> specific the file name!
1 >>> POSCAR
2 >>> CONTCAR
3 >>> CHGCAR
4 >>> LOCPOT
5 >>> maptool.cif
6 >>> maptool.xsf
7 >>> maptool.nc
8 >>> maptool.yaml
9 >>> maptool.json
--------------->>
1

    input data according to tips

     select atoms by following ways:
     1. atomic index in POSCAR
        i.e. :  1 2 4-8 10 12-30
        i.e. :  1 2 4 8 10 
     2. atomic label
        i.e. :  Si  O
     3. atomic position
        i.e. :  0 0.5 | 0.2 0.4 | 0.3 0.7
        this means atoms with 0<x<0.5, 
        0.2<y<0.4 and 0.3<z<0.7 will be seleted
        or just specific the z coordinates,
        i.e. :  ||0.3 0.7
        
--------------->>
your choice ?
--------------->>
a6
your choice ?
1 >>> build rippled structure
2 >>> build multi-layered structure
3 >>> split multi-layered structure
4 >>> resize vacuum layer
5 >>> center atomic-layer along z direction
6 >>> apply strain along different direction
7 >>> constrain atom in specific range
8 >>> get a substrate for 2D material (online!!!)
--------------->>
7
your choice ?
0 >>> specific the file name!
1 >>> POSCAR
2 >>> CONTCAR
3 >>> CHGCAR
4 >>> LOCPOT
5 >>> maptool.cif
6 >>> maptool.xsf
7 >>> maptool.nc
8 >>> maptool.yaml
9 >>> maptool.json
--------------->>
1

    input data according to tips

     select atoms by following ways:
     1. atomic index in POSCAR
        i.e. :  1 2 4-8 10 12-30
        i.e. :  1 2 4 8 10 
     2. atomic label
        i.e. :  Si  O
     3. atomic position
        i.e. :  0 0.5 | 0.2 0.4 | 0.3 0.7
        this means atoms with 0<x<0.5, 
        0.2<y<0.4 and 0.3<z<0.7 will be seleted
        or just specific the z coordinates,
        i.e. :  ||0.3 0.7
        
--------------->>
0 0.5 | 0.2 0.4 | 0.3 0.7
```
从POSCAR读取结构文件，将 ```$0<x<0.5,0.2<y<0.4,0.3<z<0.7$``` 范围内的原子固定（分数坐标）

- ### get a substrate for 2D material (online!!!)
用于寻找2D材料的可能的合成衬底,参考文献 [Penta-```$Pt_2N_4$```: an ideal two-dimensional material for nanoelectronics](http://pubs.rsc.org/en/Content/ArticleLanding/2018/NR/c8nr05561k)图i
![Pt2N4](https://github.com/haidi-ustc/image/blob/master/pmg_substrate.png?raw=true)
##### 输入文件：

POSCAR   层状材料结构文件，需要联网运行

##### 命令行参数:

```
your choice ?
--------------->>
a6
your choice ?
1 >>> build rippled structure
2 >>> build multi-layered structure
3 >>> split multi-layered structure
4 >>> resize vacuum layer
5 >>> center atomic-layer along z direction
6 >>> apply strain along different direction
7 >>> constrain atom in specific range
8 >>> get a substrate for 2D material (online!!!)
--------------->>
8
your choice ?
1 >>> input 2D structure from local disk
2 >>> get 2D structure online
--------------->>
1
your choice ?
0 >>> specific the file name!
1 >>> POSCAR
2 >>> CONTCAR
3 >>> CHGCAR
4 >>> LOCPOT
5 >>> maptool.cif
6 >>> maptool.xsf
7 >>> maptool.nc
8 >>> maptool.yaml
9 >>> maptool.json
--------------->>
1
```
从本地读取POSCAR文件，并从materialsproject数据库的81种衬底种查找与该结构匹配的衬底。

##### 输出文件:
substrate.csv 查找到的结果按照匹配度顺序保存在该文件里面，格式如下：
```
sub_id,area,film_orient,orient,sub_form
mp-1029,83.86932278192404,1 0 0,1 0 -1,BaF2
mp-1029,225.89189391462102,1 0 1,1 1 -1,BaF2
mp-1029,100.85895915643907,0 0 1,1 1 0,BaF2
mp-1029,16.80982652607318,0 0 1,1 1 1,BaF2
```

# structural analysis #
结构分析工具，获取结构对称性，单胞，惯用晶胞，XRD，指纹等信息。

## structure symmetry
分析结构的空间群信息

##### 输入文件：

random_spg_162.vasp (POSCAR  格式)

##### 命令行参数:
```
your choice ?
--------------->>
b1
your choice ?
1 >>> crystal
2 >>> molecule/cluster
--------------->>
1
your choice ?
0 >>> specific the file name!
1 >>> POSCAR
2 >>> CONTCAR
3 >>> CHGCAR
4 >>> LOCPOT
5 >>> maptool.cif
6 >>> maptool.xsf
7 >>> maptool.nc
8 >>> maptool.yaml
9 >>> maptool.json
--------------->>
0
input the full file name
--------------->>
random_spg_162.vasp
Structure Type : periodicity
Lattice Type : hexagonal
Space Group ID : 191
International Symbol : P6/mmm
Hall Symbol : -P 6 2
```

从```random_spg_162.vasp```文件读取结构文件，分析输出对称性信息：
```
Structure Type : periodicity
Lattice Type : hexagonal
Space Group ID : 191
International Symbol : P6/mmm
Hall Symbol : -P 6 2
```

## structure finger print

计算晶体或分子结构指纹信息，用于区分不同的结构，参考文献 [How to quantify energy landscapes of solids](https://aip.scitation.org/doi/10.1063/1.3079326)

###### 输入文件：
random_spg_162.vasp (POSCAR  格式)

###### 命令行参数:

```
示范法
```
## structure difference

## get primitive cell
转换结构为原胞形式

##### 输入文件：

POSCAR

##### 命令行参数:
```
your choice ?
--------------->>
b4
your choice ?
0 >>> specific the file name!
1 >>> POSCAR
2 >>> CONTCAR
3 >>> CHGCAR
4 >>> LOCPOT
5 >>> maptool.cif
6 >>> maptool.xsf
7 >>> maptool.nc
8 >>> maptool.yaml
9 >>> maptool.json
--------------->>
1
```

##### 输出文件：

maptool_primitive.vasp


## get conventional cell
将原胞转换为惯用晶胞

##### 输入文件：

POSCAR

##### 命令行参数:
```
your choice ?
--------------->>
b5
your choice ?
0 >>> specific the file name!
1 >>> POSCAR
2 >>> CONTCAR
3 >>> CHGCAR
4 >>> LOCPOT
5 >>> maptool.cif
6 >>> maptool.xsf
7 >>> maptool.nc
8 >>> maptool.yaml
9 >>> maptool.json
--------------->>
1
```

##### 输出文件：

maptool_convention.vasp


## get XRD pattern
计算结构的X射线衍射图像

CsCl 理论模拟图像
![CsCl](https://github.com/haidi-ustc/image/blob/master/CsCl_calc_XRD.png?raw=true)
CsCl 实验测量数据图像
![CsCl](https://github.com/haidi-ustc/image/blob/master/CsCl_exp_XRD.png?raw=true)

##### 输入文件：

mp-22865.cif  （CsCl晶体结构)

##### 命令行参数:

```
your choice ?
--------------->>
b6
your choice ?
0 >>> specific the file name!
1 >>> POSCAR
2 >>> CONTCAR
3 >>> CHGCAR
4 >>> LOCPOT
5 >>> maptool.cif
6 >>> maptool.xsf
7 >>> maptool.nc
8 >>> maptool.yaml
9 >>> maptool.json
--------------->>
0
input the full file name
--------------->>
mp-22865.cif
```

##### 输出文件

XRD.dat   # ```$2\theta, I$```数据

XRD.json  # json格式数据

XRD.png   # XRD图片格式文件


# vasp in/out tools 
自动产生vasp的输入文件，分析vasp输出文件

测试POSCAR
```
 FexP      
1.0
   4.40571393886535       0.000000000000000E+000  0.000000000000000E+000
  -2.20285682350647        3.81546027711513       0.000000000000000E+000
  7.186074487006176E-007  1.244664556879151E-006   9.51824425572249     
Fe P
           2           2
Direct
  0.333299994468689       0.666700005531311       0.500000000000000     
  0.666700005531311       0.333299994468689       0.500000000000000     
  0.000000000000000E+000  0.000000000000000E+000  0.000000000000000E+000
  0.000000000000000E+000  0.000000000000000E+000  0.500000000000000 
```


## prepare input files
产生特定类型的计算输入文件，可以一次性产生INCAR，POSCAR，KPOINTS，POTCAR，也可以分开产生


###  prepare all files automatically
产生和Materialsproject 数据一致性的输入文件。例如为了方便和materialsproject 里面的数据对比相稳定性，结构需要在相同的标准下计算，此时可以选择下面的输入文件。

这里的 **MIT, MP, MVL** 分别为3个标准输入文件集对应的名称


##### 输入文件：

这里并没有明确指定输入文件，因为这里用于vasp计算，所以强制要求当前目录下存在POSCAR文件

##### 命令行参数：

```
your choice ?
--------------->>
c1
------------------------ prepare intput files ------------------------
your choce ?
1 >>> prepare all files automatically
2 >>> prepare INCAR file
3 >>> prepare KPOINTS file
4 >>> prepare POTCAR file
--------------->>
1
------------------------ generate input files ------------------------
your choice?
1  >>> MIT Relax Set
2  >>> MIT NEB Set
3  >>> MIT MD Set
4  >>> MP Relax Set
5  >>> MP HSE Relax Set
6  >>> MP None SCF Set
7  >>> MP SOC Set
8  >>> MVL Elastic Set
9  >>> MVL GW Set
10 >>> MVL Slab Set
11 >>> MVL GB Set
12 >>> MVL NPT Set
--------------->>
4
MPRelaxSet(struct)

```
这里按照MP输入集产生结构优化的输入文件

##### 输出文件
input文件夹，包括四个基本输入文件，以及一个cif格式的结构文件
```
input/
├── Fe2P2.cif
├── INCAR
├── KPOINTS
├── POSCAR
└── POTCAR
```

### prepare INCAR file
自定义INCAR文件

##### 命令行参数：
```
your choice ?
--------------->>
c1
------------------------ prepare intput files ------------------------
your choce ?
1 >>> prepare all files automatically
2 >>> prepare INCAR file
3 >>> prepare KPOINTS file
4 >>> prepare POTCAR file
--------------->>
2
+----------------------------- Warning ------------------------------+
                                                                      
                        POTCAR file not found                         
                                                                      
+--------------------------------------------------------------------+
+------------------------------- Tips -------------------------------+
                                                                      
 for every letter you can append another letters for extra parameters 
  The corresponding list are:                                         
  a:  SPIN                                                            
  b:  SOC                                                             
  c:  HSE                                                             
  d:  DIPOLE correction                                               
  e:  Electric filed                                                  
  f:  Add grid                                                        
  g:  Add Pressure                                                    
  h:  DFT-D2                                                          
  i:  DFT-D3                                                          
  j:  VDW-DF                                                          
  k:  opt-B86                                                         
  l:  opt-B88                                                         
  m:  LDA+U                                                           
                                                                      
  For exmaple: aai means one optimization by condsidering SPIN and    
  DFT-D3 correction.                                                  
                                                                      
+--------------------------------------------------------------------+
------------------------ generate INCAR file -------------------------
your choice?
a >>> Optimization calculation
b >>> SCF calculation
c >>> BAND structure calculation
d >>> DOS calculation
e >>> ELF calculation
f >>> Bader charge calculation
g >>> AIMD NPT calculation
h >>> AIMD NVT calculation
i >>> Potential calculation
j >>> Partial charge calculation
k >>> STM image calculation
l >>> optical properties calculation
m >>> Mechanical properties calculation
n >>> Frequency calculation
o >>> Transition state calculation
p >>> Phonopy + vasp DFPT calculation
q >>> Phonopy + vasp finite difference calculation
--------------->>
aaf
```

这里的输入字符```aaf```中，第一个表示进行何种计算，可选范围a-q,包括优化，自洽，DOS计算，态密度计算，电子局域函数计算，bader分析，分子动力学模拟，计算势，部分电荷，STM模拟，光学性质，力学性质，频率计算，声子计算；后面的两个字母对应上面的提示项目，表示在计算过程中是否开自旋，是否加VDW修正，是否加U，是否加SOC等等。
所以这里的abf表示优化计算，打开自旋并且增加网格。

##### 输出文件：
INCAR  # 初始磁矩按照[materialsproject.org](materialsproject.org) 参考数据产生
```
#Start parameters for this run
SYSTEM  =  maptool
NWRITE  =  1
PREC    =  Accurate
ISART   =  0
ICHARG  =  2
#Electronic Relaxation 1
ENCUT   =  750.0
NELM    =  200
NELMIN  =  6
NELMDL  =  -5
EDIFF   =  1e-05
LREAL   =  Auto
#Electronic Relaxation 2
#AMIN     = 0.1
#AMIX     = 0.4
#AMIX_MAG = 1.6
#BMIX     = 1.0
#BMIX_MAG = 1.0
ALGO  =  Normal
#Ionic relaxation
EDIFFG  =   -0.01
ISIF    =    3
IBRION  =    2
POTIM   =    0.3
ISYM    =    2
NSW     =  200
# DOS related values
#EMIN     = -20.00
#EMAX     =  20.00
ISMEAR  =  0
SIGMA   =  0.05
# Write flags
LWAVE   =  False
LCHGARG  =  False
LVTOT   =  False
LVHAR   =  False
LELF    =  False
LAECHG  =  False
# Spin related parameters
# Guess values are obtained from pymatgen MPRelaxSet
ISPIN   =  2
MAGMOM  =  2*5.0 2*0.6
# parameters for add meshgrid
ADDGRID = True
```
### prepare KPOINTS file
自定义KPOINTS文件

##### 输入文件：

这里并没有明确指定输入文件，因为这里用于vasp计算，所以强制要求当前目录下存在POSCAR文件

如下命令可以产生均匀网格，用于优化或者自洽计算：
##### 命令行参数：
```
your choice ?
--------------->>
c1
------------------------ prepare intput files ------------------------
your choce ?
1 >>> prepare all files automatically
2 >>> prepare INCAR file
3 >>> prepare KPOINTS file
4 >>> prepare POTCAR file
--------------->>
3
----------------------- generate KPOINTS file ------------------------
your choice?
1 >>> automatic k-grid 
2 >>> Band structure k-path
3 >>> HSE06 k-grid
4 >>> 3D plot k-grid
--------------->>
1
 input the dimensionality and mesh grid density 
 dimensionality can be 0D 1D 2D 3D
 500 for low grid density
 1000 for medium grid density
 2000 for high grid density
 3000 for accurate density
 input format: 1 1000
--------------->>
3 1000
```
网格密度参数格式
```
体系维度   网格密度
```

##### 输出文件：
KPOINTS
```
pymatgen 4.7.6+ generated KPOINTS with grid density = 1000 / atom
0
Gamma
8 8 3
```

如下命令可以产生线性模式K点，用于能带你计算：

> 注：maptool完全按照pymatgen内部方法产生高对称点，有些高对称点可能有误，务必注意。

##### 命令行参数：
```
your choice ?
--------------->>
c1
------------------------ prepare intput files ------------------------
your choce ?
1 >>> prepare all files automatically
2 >>> prepare INCAR file
3 >>> prepare KPOINTS file
4 >>> prepare POTCAR file
--------------->>
3
----------------------- generate KPOINTS file ------------------------
your choice?
1 >>> automatic k-grid 
2 >>> Band structure k-path
3 >>> HSE06 k-grid
4 >>> 3D plot k-grid
--------------->>
2
/home/vasp/.pyenv/versions/anaconda3-4.3.0/envs/mpt_a3/lib/python3.6/site-packages/pymatgen/symmetry/bandstructure.py:63: UserWarning: The input structure does not match the expected standard primitive! The path can be incorrect. Use at your own risk.
  warnings.warn("The input structure does not match the expected standard primitive! 
```
因为与标准单胞不匹配，所以产生了警告信息

##### 输出文件：
KPOINTS
```
Line_mode KPOINTS file
30
Line_mode
Reciprocal
0.0 0.0 0.0 ! \Gamma
0.5 0.0 0.0 ! M

0.5 0.0 0.0 ! M
0.3333333333333333 0.3333333333333333 0.0 ! K

0.3333333333333333 0.3333333333333333 0.0 ! K
0.0 0.0 0.0 ! \Gamma

0.0 0.0 0.0 ! \Gamma
0.0 0.0 0.5 ! A

0.0 0.0 0.5 ! A
0.5 0.0 0.5 ! L

0.5 0.0 0.5 ! L
0.3333333333333333 0.3333333333333333 0.5 ! H

0.3333333333333333 0.3333333333333333 0.5 ! H
0.0 0.0 0.5 ! A

0.5 0.0 0.5 ! L
0.5 0.0 0.0 ! M

0.3333333333333333 0.3333333333333333 0.0 ! K
0.3333333333333333 0.3333333333333333 0.5 ! H
```


### prepare POTCAR file
自定义POTCAR文件

> 务必在软件安装阶段设置好vasp赝势搜索路径

##### 命令行参数：
```
your choice ?
--------------->>
c1
------------------------ prepare intput files ------------------------
your choce ?
1 >>> prepare all files automatically
2 >>> prepare INCAR file
3 >>> prepare KPOINTS file
4 >>> prepare POTCAR file
--------------->>
4
+------------------------------- Tips -------------------------------+
                                                                      
                     Available Pseudo-potentials are:                 
         PBE PBE_52 PBE_54 LDA LDA_52 LDA_54 PW91 LDA_US PW91_US      
                                                                      
+--------------------------------------------------------------------+
your choice ?
--------------->>
LDA
```

##### 输出文件：
产生LDA类型的POTCAR

## analysis output files
提取vasp输出文件信息，转化为能够可视化的数据

> 以下测试以Si为例子说明,POSCAR如下
```
Si2
1.0
3.348898 0.000000 1.933487
1.116299 3.157372 1.933487
0.000000 0.000000 3.866975
Si
2
direct
0.875000 0.875000 0.875000 Si
0.125000 0.125000 0.125000 Si
```

### total density of states
从vasprun.xml文件提取总态密度

##### 命令行参数：
```
your choice ?
--------------->>
c2
------------------------ vasp output analysis ------------------------
1  >>> total density of states
2  >>> projected density of states
3  >>> band structure
4  >>> projected band structure
5  >>> select one band structure
6  >>> charge density
7  >>> spin density
8  >>> charge density difference
9  >>> spin density component: up/down
10 >>> average charge density/potential
11 >>> optics analysis
12 >>> mechanical analysis
13 >>> ab initio molecular dynamics analysis
--------------->>
1
```
提取Si计算的总DOS和积分DOS

##### 屏幕输出：
```
-->> (1) Reading Data From vasprun.xml File ...
-->>     This Is a Non-Spin Calculation.
-->> (2) Writting TDOS.dat File ...
-->> (3) Writting IDOS.dat File ...
-->> (4) Saving Plot to TotalDOS.png IntegratedDOS.png File ...
```

##### 输出文件：

IDOS.dat  # 积分态密度

TDOS.dat  # 总态密度

IntegratedDOS.png # 积分态密度图
![Si IDOS](https://github.com/haidi-ustc/image/blob/master/Si_IntegratedDOS.png?raw=true)
TotalDOS.png     # 总态密度图
![Si TDOS](https://github.com/haidi-ustc/image/blob/master/Si_TotalDOS.png?raw=true)

### projected density of states
从vasprun.xml文件提取总态密度

##### 命令行参数：
```
your choice ?
--------------->>
c2
------------------------ vasp output analysis ------------------------
1  >>> total density of states
2  >>> projected density of states
3  >>> band structure
4  >>> projected band structure
5  >>> select one band structure
6  >>> charge density
7  >>> spin density
8  >>> charge density difference
9  >>> spin density component: up/down
10 >>> average charge density/potential
11 >>> optics analysis
12 >>> mechanical analysis
13 >>> ab initio molecular dynamics analysis
--------------->>
2
-->> (1) Reading Data From vasprun.xml File ...
-->> (2) Reading Data From PROCAR File ...

    input data according to tips

     select atoms by following ways:
     1. atomic index in POSCAR
        i.e. :  1 2 4-8 10 12-30
        i.e. :  1 2 4 8 10 
     2. atomic label
        i.e. :  Si  O
     3. atomic position
        i.e. :  0 0.5 | 0.2 0.4 | 0.3 0.7
        this means atoms with 0<x<0.5, 
        0.2<y<0.4 and 0.3<z<0.7 will be seleted
        or just specific the z coordinates,
        i.e. :  ||0.3 0.7
        
--------------->>
1
-->>     This Is a Non-Spin Calculation.
-->> (3) Writting Projected DOS Data to PDOS.dat File ..
```
计算第一个原子对于DOS的贡献


##### 屏幕输出：
```
-->>     This Is a Non-Spin Calculation.
-->> (3) Writting Projected DOS Data to PDOS.dat File ...
```

##### 输出文件：
PDOS.dat,文件格式

```
#String: 1
#Selected atom: 1
#  K-Distance  Energy(ev)           s          py          pz          px         dxy         dyz         dz2         dxz         dx2
  -13.144153    0.000000    0.000000    0.000000    0.000000    0.000000    0.000000    0.000000    0.000000    0.000000
  -13.066453    0.000000    0.000000    0.000000    0.000000    0.000000    0.000000    0.000000    0.000000    0.000000
  -12.988653    0.000000    0.000000    0.000000    0.000000    0.000000    0.000000    0.000000    0.000000    0.000000
  -12.910853    0.000000    0.000000    0.000000    0.000000    0.000000    0.000000    0.000000    0.000000    0.000000
  -12.833053    0.000000    0.000000    0.000000    0.000000    0.000000    0.000000    0.000000    0.000000    0.000000
  -12.755353    0.000000    0.000000    0.000000    0.000000    0.000000    0.000000    0.000000    0.000000    0.000000
  -12.677553    0.000000    0.000000    0.000000    0.000000    0.000000    0.000000    0.000000    0.000000    0.000000
```

### band structure
从vasprun.xml读取能带数据，同时从KPOINTS读取标签信息，从OUTCAR读取磁矩信息

##### 命令行参数：
```
your choice ?
--------------->>
c2
------------------------ vasp output analysis ------------------------
1  >>> total density of states
2  >>> projected density of states
3  >>> band structure
4  >>> projected band structure
5  >>> select one band structure
6  >>> charge density
7  >>> spin density
8  >>> charge density difference
9  >>> spin density component: up/down
10 >>> average charge density/potential
11 >>> optics analysis
12 >>> mechanical analysis
13 >>> ab initio molecular dynamics analysis
--------------->>
3
```
##### 屏幕输出：
```
-->> (1) Reading Data From vasprun.xml File ...
-->> (2) Reading Data From KPOINTS File ...
-->> (3) Reading Data From OUTCAR File ...
-->>     This Is a Non-Spin Calculation.
-->>     This Material Is a Semiconductor.
-->>     vbm=5.616700 eV cbm=6.227100 eV gap=0.610400 eV
-->> (4) Writting Band Structure Data to BAND.dat File ...
-->> (5) Writting Label infomation to HighSymmetricPoints.dat File ...
-->> (6) Saving Plot to BAND.png File ...
```

##### 输出文件：
BAND.dat   #能带

HighSymmetricPoints.dat  #高对称点信息
```
#       index       label    position
           0    $\Gamma$    0.000000
           1           X    1.148930
           2           X    1.148930
           3           W    1.723395
           4           W    1.723395
           5           K    2.129603
           6           K    2.129603
           7    $\Gamma$    3.348227
           8    $\Gamma$    3.348227
           9           L    4.343229
          10           L    4.343229
          11           U    5.046802
          12           U    5.046802
          13           W    5.453010
          14           W    5.453010
          15           L    6.265426
          16           L    6.265426
          17    K$\mid$U    6.968999
          18           X    7.375207
```

BAND.png #能带图
![Si band](https://github.com/haidi-ustc/image/blob/master/Si_band.png?raw=true)


### projected band structure
从vasprun.xml, PROCAR以及KPIOINTS文件提取投影能带信息

##### 命令行参数：
```
your choice ?
--------------->>
c2
------------------------ vasp output analysis ------------------------
1  >>> total density of states
2  >>> projected density of states
3  >>> band structure
4  >>> projected band structure
5  >>> select one band structure
6  >>> charge density
7  >>> spin density
8  >>> charge density difference
9  >>> spin density component: up/down
10 >>> average charge density/potential
11 >>> optics analysis
12 >>> mechanical analysis
13 >>> ab initio molecular dynamics analysis
--------------->>
4
-->> (1) Reading Data From vasprun.xml File ...
-->> (2) Reading Data From PROCAR File ...
-->> (3) Reading Data From KPOINTS File ...

    input data according to tips

     select atoms by following ways:
     1. atomic index in POSCAR
        i.e. :  1 2 4-8 10 12-30
        i.e. :  1 2 4 8 10 
     2. atomic label
        i.e. :  Si  O
     3. atomic position
        i.e. :  0 0.5 | 0.2 0.4 | 0.3 0.7
        this means atoms with 0<x<0.5, 
        0.2<y<0.4 and 0.3<z<0.7 will be seleted
        or just specific the z coordinates,
        i.e. :  ||0.3 0.7
        
--------------->>
2
```
这里选择第二个原子的投影贡献


##### 屏幕输出：
```
-->>     This Is a Non-Spin Calculation.
-->> (4) Writting Projected Band Structure Data to PBAND.dat File ...
-->> (5) Writting Label infomation to HighSymmetricPoints.dat File ...
```
##### 输出文件：
PBAND.dat #投影数据，这里设置LORBIT=10，故只投影到spd轨道，格式如下
```
#String: 1
#Selected atom: 2
#  K-Distance  Energy(ev)           s           p           d
    0.000000  -12.085035    0.212000    0.000000    0.000000
    0.039618  -12.079935    0.000000    0.214000    0.000000
    0.079237  -12.064535    0.000000    0.214000    0.000000
    0.118855  -12.039035    0.000000    0.214000    0.000000
    0.158473  -12.003235    0.000000    0.151000    0.000000
    0.198091  -11.957335    0.000000    0.151000    0.000000
    0.237710  -11.901235    0.000000    0.151000    0.000000
```
 HighSymmetricPoints.dat  #高对称点信息，和band计算一样

### select one band structure
选择特定的能带，用于对此能带画图是进行强调分析或者用此条能带进行有效质量的计算

##### 命令行参数：
```
your choice ?
--------------->>
------------------------ vasp output analysis ------------------------
1  >>> total density of states
2  >>> projected density of states
3  >>> band structure
4  >>> projected band structure
5  >>> select one band structure
6  >>> charge density
7  >>> spin density
8  >>> charge density difference
9  >>> spin density component: up/down
10 >>> average charge density/potential
11 >>> optics analysis
12 >>> mechanical analysis
13 >>> ab initio molecular dynamics analysis
--------------->>
5
-->> (1) Reading Data From vasprun.xml File ...
-->> (2) Reading Data From KPOINTS File ...
-->>     This Is a Non-Spin Calculation.
-->>     Total band number is 8
-->>     Total electron number is 8.0
which band would like to select ?
--------------->>
4
```
命令行提示Si体系具有8条能带，总电子数为8，这里选择第四条能带，即价带

##### 输出文件：

BAND_4.dat 


### charge density
提取总电荷密度文件(无论自旋极化与否)，VESTA可视化

##### 命令行参数，
```
your choice ?
--------------->>
c2
------------------------ vasp output analysis ------------------------
1  >>> total density of states
2  >>> projected density of states
3  >>> band structure
4  >>> projected band structure
5  >>> select one band structure
6  >>> charge density
7  >>> spin density
8  >>> charge density difference
9  >>> spin density component: up/down
10 >>> average charge density/potential
11 >>> optics analysis
12 >>> mechanical analysis
13 >>> ab initio molecular dynamics analysis
---------
6
```

##### 输出文件：
CHARGE.vasp


### spin density
自旋密度：上自旋-下自旋，得到净电荷密度，用于观察磁矩分布，VESTA可视化
```
your choice ?
--------------->>
c2
------------------------ vasp output analysis ------------------------
1  >>> total density of states
2  >>> projected density of states
3  >>> band structure
4  >>> projected band structure
5  >>> select one band structure
6  >>> charge density
7  >>> spin density
8  >>> charge density difference
9  >>> spin density component: up/down
10 >>> average charge density/potential
11 >>> optics analysis
12 >>> mechanical analysis
13 >>> ab initio molecular dynamics analysis
---------
7
```

#####  输出文件：
Spin.vasp  # 自旋电荷密度

Total.vasp  # 总电荷密度


### charge density difference
计算电荷密度之差，主要用于计算分子吸附在固体表面之后导致的电荷转移情况(定性，VESTA可视化)

##### 输入文件：
两个或多个电荷密度文件

##### 命令行参数：

```
your choice ?
--------------->>
c2
------------------------ vasp output analysis ------------------------
1  >>> total density of states
2  >>> projected density of states
3  >>> band structure
4  >>> projected band structure
5  >>> select one band structure
6  >>> charge density
7  >>> spin density
8  >>> charge density difference
9  >>> spin density component: up/down
10 >>> average charge density/potential
11 >>> optics analysis
12 >>> mechanical analysis
13 >>> ab initio molecular dynamics analysis
--------------->>
8
input the name of charge density file like this
file_A file_B, or
it means rho(A)-rho(B))
file_A file_B file_C
it means rho(A)-rho(B)-rho(C)
--------------->>
CHG CHG CHG
```

这里```CHG CHG CHG```表示 ```CHG-CHG-CHG=-CHG```仅仅用来测试，在实际计算过程中可能需要```$\rho(subsrate+mol)-\rho(subsrate)-\rho(mol)$```,只需要顺次输入相应的文件名即可即可。

##### 屏幕输出：
```
-->> (1) Reading Charge Density From CHG File ...
-->> (2) Reading Charge Density From CHG File ...
-->> (3) Reading Charge Density From CHG File ...
-->> (4) Writing Charge Density Difference to Diff.vasp  File ...
```
##### 输出文件：
Diff.vasp 


### spin density component: up/down
对自旋极化计算提取自旋向上和自旋向下的电荷密度

##### 命令行参数：
```
your choice ?
--------------->>
c2
------------------------ vasp output analysis ------------------------
1  >>> total density of states
2  >>> projected density of states
3  >>> band structure
4  >>> projected band structure
5  >>> select one band structure
6  >>> charge density
7  >>> spin density
8  >>> charge density difference
9  >>> spin density component: up/down
10 >>> average charge density/potential
11 >>> optics analysis
12 >>> mechanical analysis
13 >>> ab initio molecular dynamics analysis
--------------->>
9
```
##### 命令行输出：
```
-->> (1) Reading Spin Up/Down Density From CHG File ...
-->> (2) Writing Spin Up Density To Spinup.vasp File ...
-->> (3) Writing Spin Down Density To Spindown.vasp File ..
```
##### 输出文件：
Spinup.vasp     #自旋向下

spindown.vasp  # 自旋向下

### average charge density/potential

平均电荷密度或者势函数

##### 命令行输入：
```
your choice ?
--------------->>
c2
------------------------ vasp output analysis ------------------------
1  >>> total density of states
2  >>> projected density of states
3  >>> band structure
4  >>> projected band structure
5  >>> select one band structure
6  >>> charge density
7  >>> spin density
8  >>> charge density difference
9  >>> spin density component: up/down
10 >>> average charge density/potential
11 >>> optics analysis
12 >>> mechanical analysis
13 >>> ab initio molecular dynamics analysis
--------------->>
10
which file would like to average: CHG or LOCPOT ?
--------------->>
CHG
-->> (1) Reading Data From CHG File ...
which direction would like to average: x y or z ?
--------------->>
x
```
沿着x方向对电荷密度进行平均

##### 输出文件：
average_CHG_x.dat



### optics analysis
线性光学性质计算

##### 命令行参数：
```
your choice ?
--------------->>
c2
------------------------ vasp output analysis ------------------------
1  >>> total density of states
2  >>> projected density of states
3  >>> band structure
4  >>> projected band structure
5  >>> select one band structure
6  >>> charge density
7  >>> spin density
8  >>> charge density difference
9  >>> spin density component: up/down
10 >>> average charge density/potential
11 >>> optics analysis
12 >>> mechanical analysis
13 >>> ab initio molecular dynamics analysis
--------------->>
11
```
##### 输出文件：
AbsorbSpectrum.dat  

EnergyLossSpectrum.dat  

ExtictionSpectrum.dat  

ReflectivitySpectrum.dat  

RefractiveSpectrum.dat

### mechanical analysis

### ab initio molecular dynamics analysis

暂不支持


#  online retrieve 
在线从[materialsproject.org](https://www.materialsproject.org)查找数据库下载数据，包括结构，能带，态密度，基本的元素信息，能量信息等。

## get band/dos by mp-ID
通过mp-ID查找数据

如图为mp-3748的能带结构图

![mp-3748](https://github.com/haidi-ustc/image/blob/master/band.png?raw=true)

##### 命令行参数：

```
your choice ?
--------------->>
e1
input the mp-ID
--------------->>
mp-123
```
用于下载id号为mp-123的数据

##### 屏幕输出：

```
ComputedEntry mp-123 - Nd4
Energy = -19.0604
Correction = 0.0000
Parameters:
run_type = GGA
is_hubbard = False
pseudo_potential = {'functional': 'PBE', 'labels': ['Nd_3'], 'pot_type': 'paw'}
hubbards = {}
potcar_symbols = ['PBE Nd_3']
oxide_type = None
Data:
oxide_type = None
```



##### 输出文件:

mp-123_band.png    #能带图
![band](https://github.com/haidi-ustc/image/blob/master/mp-123_band.png?raw=true)
mp-123_dos.png     #态密度图
![band](https://github.com/haidi-ustc/image/blob/master/mp-123_dos.png?raw=true)


## get structure from materialsproject database
按照特定的查询标准下载结构文件，提供三种查询方法：
- mp-ID 材料ID编号
```
mp-123
```
- 材料的分子式 
```
Fe3O4
```
或者
```
(SiO2)2
```
- 材料所含元素 （用空格分开即可）
```
Si O Ca
```

##### 命令行参数:
```
your choice ?
--------------->>
e2
your choice ?
1 >>> get a structure by mp-ID
2 >>> get a structure by fomular
3 >>> get a structure by elements
--------------->>
1
input the mp-ID
--------------->>
mp-123
```
指定下载mp-ID为mp-123的结构

##### 输出文件：

mp-123.vasp



## get properties by mp-ID
下载指定ID对应的所有数据信息


##### 命令行参数：

```
your choice ?
--------------->>
e3
input the mp-ID
--------------->>
mp-123
```

##### 输出文件：

mp-123.json  # json 格式数据，包含结构的能量，带隙，对称性等信息。内容如下：
```
[
    {
        "energy": -19.06038229,
        "energy_per_atom": -4.7650955725,
        "volume": 141.75484566720664,
        "formation_energy_per_atom": 0.0,
        "nsites": 4,
        "unit_cell_formula": {
            "Nd": 4.0
        },
        "pretty_formula": "Nd",
        "is_hubbard": false,
        "elements": [
            "Nd"
        ],
        "nelements": 1,
        "e_above_hull": 0,
        "hubbards": {},
        "is_compatible": true,
        "spacegroup": {
            "source": "spglib",
            "symbol": "P6_3/mmc",
            "number": 194,
            "point_group": "6/mmm",
            "crystal_system": "hexagonal",
            "hall": "-P 6c 2c"
        },
        "task_ids": [
            "mp-123",
            "mp-905471",
            "mp-918603",
            "mp-919579"
        ],
        "band_gap": 0.0,
        "density": 6.758695862009324,
        "icsd_id": null,
        "icsd_ids": [
            645581,
            645585,
            76592,
            102657,
            645584,
            43577,
            164281,
            645579,
            645577,
            150953
        ],
        "total_magnetization": 0.0015692,
        "material_id": "mp-123",
        "oxide_type": "None",
        "tags": [
            "High pressure experimental phase",
            "Neodymium - alpha",
            "Neodymium"
        ],
        "elasticity": null,
        "full_formula": "Nd4"
    }
]
```
> 上述输出文件中cif部分由于内容过多，这没有显示






