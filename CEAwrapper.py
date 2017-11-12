# -*- coding: utf-8 -*-
import os
import sys
import json
import numpy as np
from scipy import interpolate
import subprocess
import warnings
# warnings.filterwarnings('ignore')

class Fuel:
    def __init__(self, name, mass_ratio, enthalpy, num_C, num_H, num_O):
        self.name = name
        self.mass_ratio = mass_ratio
        self.enthalpy = enthalpy
        self.num_C = num_C
        self.num_H = num_H
        self.num_O = num_O


class CEAwrapper:
    def __init__(self, cea_config_file, AR, Pair = 0.1013):
        self.initialize()
        self.fuel_read(cea_config_file)
        self.AR = AR
        self.Pair = Pair

    def runCEA(self, input_file):
        if os.name == 'nt':
            cmd = './FCEA2.exe'
        else:
            cmd = './FCEA2'
        p = subprocess.Popen(cmd, shell=False, stdin=subprocess.PIPE, stdout=subprocess.PIPE)
        # さっき作った.inp ファイルの，".inp"を除いたファイル名をCEAの標準入力に渡す
        # なぜ argv で渡せるようにしてくれなかったのか
        if os.name == 'nt':
            p.communicate(input_file + b'\n') # Win環境ではコマンド末尾に¥nを付けないとバグる
        else:
            p.communicate(input_file)
        p.wait()

    def initialize(self):
        os.chdir('./CEA')
        if os.name == 'nt':
            if os.path.isfile('FCEA2.exe') and os.path.isfile('b1b2b3.exe') and os.path.isfile('syntax.exe'):
                pass
            else:
                p = subprocess.Popen('make', shell=False, stdin=subprocess.PIPE, stdout=subprocess.PIPE)
                p.wait()
        else:
            if os.path.isfile('FCEA2') and os.path.isfile('b1b2b3') and os.path.isfile('syntax'):
                pass
            else:
                p = subprocess.Popen('make', shell=False, stdin=subprocess.PIPE, stdout=subprocess.PIPE)
                p.wait()
        if os.path.isfile('thermo.lib'):
            pass
        else:
            self.runCEA(b'thermo')
        if os.path.isfile('trans.lib'):
            pass
        else:
            self.runCEA(b'trans')
        os.chdir('../')

    def fuel_read(self, cea_config_file):
        json_obj = json.load(open(cea_config_file))
        self.fuel_list = []
        for i in range(100):
            try:
                component = json_obj.get('FuelCEA').get('Component-'+str(i))
                if json_obj.get('FuelCEA').get('Component-'+str(i)).get('exist') == True:
                    name = json_obj.get('FuelCEA').get('Component-'+str(i)).get('Name')
                    mass_ratio = json_obj.get('FuelCEA').get('Component-'+str(i)).get('Mass Ratio [%]')
                    enthalpy = json_obj.get('FuelCEA').get('Component-'+str(i)).get('Enthalpy [kJ/mol]')
                    num_C = json_obj.get('FuelCEA').get('Component-'+str(i)).get('C')
                    num_H = json_obj.get('FuelCEA').get('Component-'+str(i)).get('H')
                    num_O = json_obj.get('FuelCEA').get('Component-'+str(i)).get('O')
                    fuel = Fuel(name, mass_ratio, enthalpy, num_C, num_H, num_O)
                    self.fuel_list.append(fuel)
            except (KeyError, AttributeError):
                pass

    def cea_iac(self, OF, Pc):
        os.chdir('./CEA')
        file_name = b'ceatemp' # 一時的に作られるファイル
        # ---- .inpファイル作り ----
        input_file = file_name + b'.inp'
        f = open(input_file, 'w') # 書き込みモードで開く

        f.write('problem\n')
        f.write('o/f='+str(OF)+',\n')
        f.write('rocket frozen nfz=2\n')
        f.write('p,bar='+str(Pc*10)+'\n')
        f.write('sup,ae/at='+str(self.AR)+'\n')
        f.write('react\n')
        f.write('oxid=N2O wt=100 t,k=298.15\n')
        for fuel in self.fuel_list:
            f.write('fuel=' + str(fuel.name) + ' wt=' + str(fuel.mass_ratio) + ' t,k=298.15\n')
            f.write('h,kj/mol=' + str(fuel.enthalpy) + ' C ' + str(fuel.num_C) + ' H ' + str(fuel.num_H) + ' O ' + str(fuel.num_O) + '\n')
        f.write('output\n')
        f.write('plot p ispfz ivacfz cffz o/f son machfz\n')
        f.write('end\n')

        f.close() # ファイルを閉じる
        # plotファイルから数字を読み出すがispやcfなどは末尾にfzを付けないとfrozen条件の値ではなく
        # equilibrium 条件の値が出力されてしまうので注意．
        # iac オプションを付けると出力は inj, throat, exit の順になる．
        # exit の指定は出口圧力比(pi/p)，開口比(suparまたはsup,ae/at)などの方法があり，複数指定すると
        # exit の出力も複数になるが，ここではある決まった(作った)ノズルに対する性能を見たいので開口比で指定．

        self.runCEA(file_name)

        # ---- .pltファイル読み取り ----
        output_file = file_name + b'.plt'
        f = open(output_file, 'r')
        lines = f.readlines()
        line = lines[3].split()
        Pe = float(line[0])/10
        isp_opt = float(line[1])
        ivac = float(line[2])
        cf_opt = float(line[3])
        Pc = float(lines[1].split()[0])/10
        OF = float(line[4])
        Cs = float(line[5])
        Mach = float(line[6])
        # CEA は開口比を与えた時には大気圧によらず最適膨張(= 圧力推力なし)を仮定して
        # Ispやcfを出してくる．従って，「ある大気圧下での性能(海面高度含む)」を知りたい時は
        # 圧力推力項を加減してやる必要がある．それが以下の3式．教科書通り．
        try:
            cf = cf_opt + self.AR / Pc * (Pe - self.Pair)
            cf = max(0.1, cf)
            isp = (isp_opt * cf / cf_opt) / 9.80665
            cstar = isp * 9.80665 / cf
            Ve = Cs * Mach
        except ZeroDivisionError:
            cf = 0.1
            isp = 0.1
            cstar = 0.1
            Ve = 0.0

        # ファイル削除
        # os.remove(input_file)
        # os.remove(output_file)

        output_data = [cf, cstar, isp, Pe, Ve, Pc, OF, self.AR]
        os.chdir('../')
        return output_data


if __name__ == '__main__':
    import matplotlib.pyplot as plt
    
    AR = 3.49
    Pair = 0.1013
    
    cea = CEAwrapper('CEAconfig.json', 3.49)
    OF_array = np.arange(1.0, 11, 1)
    Pc_array = np.arange(0.3, 10.3, 1)
    OF_temp, Pc_temp = np.meshgrid(OF_array, Pc_array)
    # CF_array = list(map(cea.cea_iac, OF_temp, Pc_temp, AR, Pair))
    # cstar_array = list(map(cea.cea_iac, OF_temp, Pc_temp, AR, Pair))

    def dd(OF, Pc):
        return cea.cea_iac(OF, Pc)[1]

    data = list(map(dd, OF_temp.ravel(), Pc_temp.ravel()))

    inter = interpolate.interp2d(OF_array, Pc_array, data, kind='cubic')

    OFnew = np.array([6.3, 6.2, 7.0])
    Pcnew = np.array([2.3, 2.2, 3.0])

    print(inter(OFnew.ravel(), Pcnew.ravel()))
