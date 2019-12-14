import sqlite3

Dbpath = "C:/Users/double/Desktop/毕业设计"
Calc = {}

# test case
SmilesCase = "CC(=O)C1(O)CC2OC12";


# SmilesCase = "CCC1COC1(C)CC"

class Graph:
    def __init__(self):
        self.AllAtoms = []
        self.AllBonds = []


'''
type:
1 单键
2 双键
3 三键
'''


class Bond:
    def __init__(self, type1, index1):
        self.Type = type1
        self.index = index1
        self.Atom1 = None
        self.Atom2 = None


'''
Type: 
C 1
O 2
H 3

'''


class Atom:
    def __init__(self):
        self.Type = -1
        self.index = -1
        self.bond = []
        self.neighborAtom = []

    def __init__(self, type1, index1):
        self.Type = type1;
        self.index = index1;
        self.bond = []
        self.neighborAtom = []
        self.branchIndex = -1


G = Graph()
E = []
'''
基团类型：

1、伯仲叔季碳
C/C/H3  :11
C/C2/H2 :12
C/C3/H  :13
C/C4    :14
2、含氧基团
C/
'''


def Smiles2Vector(Smiles: str) -> []:
    # 新建一个虚的原子，作为所有原子的开始
    AtomIndex = 1
    StartAtom = Atom(-1, AtomIndex)
    StartAtom.index = 0
    preAtom = StartAtom
    G.AllAtoms.append(StartAtom)
    BranchList = {1: -1, 2: -1, 3: -1, 4: -1, 5: -1, 6: -1, 7: -1, 8: -1, 9: -1}  # 记录分支情况
    length = len(Smiles)
    if length == 0:
        return []
    branchAtom = Atom(-1, -1)
    BondIndex = 1
    # 此处为关键处理部分
    preBond = 1
    for atom in Smiles:
        if atom == 'C':
            Catom = Atom(1, AtomIndex)  # 新建原子
            AtomIndex += 1  # 原子编号更新
            preAtom.neighborAtom.append(Catom)  # 将该原子与上一个原子相连
            Catom.neighborAtom.append(preAtom)
            TempBond = Bond(preBond, BondIndex)  # 新建键
            preBond = 1  # 更新键之后默认为单键，除非读到其他键
            TempBond.Atom1 = preAtom  # 更新键两边的原子
            TempBond.Atom2 = Catom
            BondIndex += 1
            preAtom.bond.append(TempBond)
            Catom.bond.append(TempBond)
            preAtom = Catom

            G.AllAtoms.append(Catom)
            G.AllBonds.append(TempBond)
            del TempBond
            del Catom
        if atom == 'O':
            Oatom = Atom(2, AtomIndex)  # 新建原子
            AtomIndex += 1  # 原子编号更新
            preAtom.neighborAtom.append(Oatom)  # 将该原子与上一个原子相连
            Oatom.neighborAtom.append(preAtom)  # 将上一个原子与该原子相连
            TempBond = Bond(preBond, BondIndex)  # 新建键
            preBond = 1  # 更新键之后默认为单键，除非读到其他键
            TempBond.Atom1 = preAtom  # 更新键两边的原子
            TempBond.Atom2 = Oatom
            BondIndex += 1
            preAtom.bond.append(TempBond)  # 为前一个原子添加该键
            Oatom.bond.append(TempBond)  # 为当前原子添加该键
            preAtom = Oatom

            G.AllAtoms.append(Oatom)
            G.AllBonds.append(TempBond)
            del TempBond
            del Oatom
        # 此时读到一个双键
        if atom == '=':
            preBond = 2
        # 此时读到一个三键
        if atom == '#':
            preBond = 3
        # 此时读到分支开始
        if atom == '(':
            branchAtom = preAtom  # 读到分支，则记录之前的原子作为分支开始

        # 此时读到分支结束
        if atom == ')':
            preAtom = branchAtom  # 读到分支结束，则下一个原子应该与分支前的原子相连
        # 此时读到一个环
        if 48 < ord(atom) < 58:
            RingIndex = ord(atom) - 48
            if BranchList[RingIndex] == -1:
                BranchList[RingIndex] = preAtom.index  # 之前为出现成环符号，则把当前原子标记成成环符号

            else:
                PreRingIndex = BranchList[RingIndex]
                PreRingAtom = G.AllAtoms[PreRingIndex]
                assert PreRingAtom.index == PreRingIndex
                PreRingAtom.neighborAtom.append(preAtom)

                preAtom.neighborAtom.append(PreRingAtom)
                TempBond = Bond(1, BondIndex)
                BondIndex += 1
                TempBond.Atom1 = preAtom
                TempBond.Atom2 = PreRingAtom
                preAtom.bond.append(TempBond)
                PreRingAtom.bond.append(TempBond)

                G.AllBonds.append(TempBond)
                del TempBond

    # 通过Floyd算法找到点到所有点最短距离

    length = len(G.AllAtoms)
    for i in range(length):
        E.append([])
        for j in range(length):
            E[i].append(9999)

    for i in range(length):
        E[i][i] = 0
    for Atoms in G.AllAtoms:
        tempx = Atoms.index
        for NeighborAtom in Atoms.neighborAtom:
            tempy = NeighborAtom.index
            E[tempx][tempy] = 1
    # 已经构建了Floyd矩阵
    for i in range(length):
        for j in range(length):
            for k in range(length):
                if E[i][j] > E[i][k] + E[k][j]:
                    E[i][j] = E[i][k] + E[k][j]

    return []


# def BuildFloyd():


def Analysis_OCH_pre():
    # 统计数据库中的元素都有哪些
    conn = sqlite3.connect(Dbpath + '/g4mp2-gdb9.db')
    cursor = conn.cursor()

    sql = "SELECT value FROM text_key_values WHERE key = \'Smiles\'"
    cursor.execute(sql)
    Value = cursor.fetchall()
    for molecular in Value:
        for c in molecular[0]:
            if c in Calc:
                Calc[c] += 1
            else:
                Calc[c] = 1
    # 统计显示，数据库分子中只含CHONF
    # 对数据库做一个基本统计，查看数据库分子构成，挑选含CHO的条目

    return


def Analysis_OCH():
    conn = sqlite3.connect(Dbpath + '/g4mp2-gdb9.db')
    cursor = conn.cursor()

    sql = "SELECT value FROM text_key_values WHERE key = \'Smiles\'"
    cursor.execute(sql)
    Value = cursor.fetchall()
    count = 0
    All_sum = 0
    for molecule in Value:
        if 'N' in molecule[0] or 'n' in molecule[0] or 'F' in molecule[0]:
            All_sum += 1
        else:
            count += 1
            All_sum += 1
    print(count)
    print(All_sum)


def Group_Statistics():
    conn = sqlite3.connect(Dbpath + '/g4mp2-gdb9.db')
    cursor = conn.cursor()
    sql = "SELECT value FROM text_key_values WHERE key = \'Plain_Smiles\'"
    cursor.execute(sql)
    Value = cursor.fetchall()
    Group_Type = []
    global G
    global E
    Type = 0
    for Smiles in Value:
        if 'N' in Smiles[0] or 'n' in Smiles[0] or 'F' in Smiles[0]:
            pass
        else:
            Smiles2Vector(Smiles[0])
            print(Type)
            Type+=1
            for atom in G.AllAtoms:
                for bond in atom.bond:
                    AnotherAtom = None
                    if bond.Atom2.Type == atom.Type:
                        AnotherAtom = bond.Atom1
                    else:
                        AnotherAtom = bond.Atom2
            G = Graph()
            E = []
    return


def mainProcess():
    conn = sqlite3.connect(Dbpath + '/g4mp2-gdb9.db')
    cursor = conn.cursor()
    sql = "SELECT value FROM text_key_values WHERE key = \'Plain_Smiles\'"
    cursor.execute(sql)
    Value = cursor.fetchall()


'''
input is Smiles expression, output is the vector used for ML model
Smiles only contains C H O and bond structure information
Output vector 

分析思路：
读取表达式，提取基团，以树状结构组织基团，计算相关的参数
'''

if __name__ == '__main__':
    Group_Statistics()
