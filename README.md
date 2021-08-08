# CellCommunication
I imitated CellPhoneDB to achieve the communication on my data. 

# 细胞通讯  --  详见 [细胞通讯实战演练](https://www.jianshu.com/p/627d0a0dae55)
## 背景
细胞-细胞间的交流是肿瘤细胞免疫微环境一个重要方向，比如当下比较火热的免疫治疗（antiPD1 or antiPD-L1）。细胞-细胞的交流指的是，细胞A表达ligand, 细胞B表达receptor, 细胞A释放的，作用于细胞B，由此引发了细胞B一系列的生物学变化，如细胞迁移，凋亡，增殖分化等。通过实验的方式筛选**有意义的受体-配体**，人的时间和精力有限。因此我们希望能够通过，生物信息的方法，快速筛选出具有统计意义的受体-配体

## 筛选方法
- 参考文章methods: Cell-cell interaction analysis

- [CellPhoneDB官方文档](https://links.jianshu.com/go?to=https%3A%2F%2Fwww.cellphonedb.org%2Fexplore-sc-rna-seq)

## 数据类型
所用到的数据在output/下，matrix.txt是稀疏性矩阵，需要你用**R包Matrix**读入，读入之后可转化为你常见的矩阵，该矩阵行代表基因，列代表每个单细胞；mdsc_tumor_nk_celltype.txt是对细胞类型的注释。（含有3种细胞类型：分别是T细胞，MDSC细胞（骨髓来源抑制细胞），Tumor细胞）；zhang_tables7.xlsx是候补的受体-配体文件，已被划分成L-Rpairs.txt文件和geneOfpairs.txt；symbol_id.txt是基因别名的对应关系；barcodes.txt和features.txt分别代表细胞和基因

## 实现要求
- 根据配体-受体的p-value, 挑出Tumor-NK，Tumor-MDSC，NK-Tumor，NK-MDSC，MDSC-NK，MDSC-Tumor这六种细胞-细胞交流，每组最显著的5个配对。（NK-Tumor：代表T细胞释放配体，Tumor细胞表达受体，其余类似）

- 参照文章fig7.E,画出CD274-PDCD1在我们数据中的interaction强弱（NK，Tumor，MDSC）--**和弦图**
