import sqlite3

Dbpath = "C:/Users/double/Desktop/毕业设计"
Calc = {}

# test case
#SmilesCase = "CC(=O)C1(O)CC2OC12";
SmilesCase = "CC=CCO"

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
    def __init__(self,type1,index1):
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
        self.type = type1;
        self.index = index1;
        self.bond = []
        self.neighborAtom = []


G = Graph()

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
    StartAtom = Atom(-1,AtomIndex)
    StartAtom.index = 0
    preAtom = StartAtom
    G.AllAtoms.append(StartAtom)

    length = len(Smiles)
    if length == 0:
        return []
    branchAtoms = Atom(-1,-1)
    BondIndex = 1
    # 此处为关键处理部分
    preBond = 1
    for atom in Smiles:
        if atom == 'C':
            Catom = Atom(1, AtomIndex)#新建原子
            AtomIndex += 1#原子编号更新
            preAtom.neighborAtom.append(Catom)#将该原子与上一个原子相连
            Catom.neighborAtom.append(preAtom)
            TempBond = Bond(preBond,BondIndex)#新建键
            preBond = 1#更新键之后默认为单键，除非读到其他键
            TempBond.Atom1 = preAtom#更新键两边的原子
            TempBond.Atom2 = 1
            BondIndex += 1
            preAtom.bond.append(TempBond)
            Catom.bond.append(TempBond)
            preAtom = Catom

            G.AllAtoms.append(Catom)
            G.AllBonds.append(TempBond)
            del TempBond
            del Catom
        if atom == 'O':
            Oatom = Atom(2, AtomIndex)#新建原子
            AtomIndex += 1#原子编号更新
            preAtom.neighborAtom.append(Oatom)#将该原子与上一个原子相连
            Oatom.neighborAtom.append(preAtom)#将上一个原子与该原子相连
            TempBond = Bond(preBond,BondIndex)#新建键
            preBond = 1#更新键之后默认为单键，除非读到其他键
            TempBond.Atom1 = preAtom#更新键两边的原子
            TempBond.Atom2 = 2
            BondIndex += 1
            preAtom.bond.append(TempBond)#为前一个原子添加该键
            Oatom.bond.append(TempBond)#为当前原子添加该键
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
            pass
        # 此时读到分支结束
        if atom == ')':
            pass
        # 此时读到一个环
        if 0 < ord(atom) < 9:
            pass




    return []


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


'''
input is Smiles expression, output is the vector used for ML model
Smiles only contains C H O and bond structure information
Output vector 

分析思路：
读取表达式，提取基团，以树状结构组织基团，计算相关的参数
'''

if __name__ == '__main__':
    Smiles2Vector(SmilesCase)
