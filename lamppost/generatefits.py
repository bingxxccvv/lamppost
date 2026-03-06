import os
import numpy as np
from astropy.io import fits
#1.1是test
#1.2是正式1.1-100
# --- 1. 参数与路径配置 ---
folder = "/data/home/bingbing/mywork/lamppost/data/"
output_file = 'lp_thickness_v1.2.fits'

# Q 列表 (根据你的需求，这里可以遍历 Q)
Q_list = [-0.25, 0.0, 0.25, 0.375, 0.5, 0.5625, 0.625, 0.6875, 0.75, 0.8125, 0.875, 0.9375]

# 归一化系数 (确保长度能覆盖你的 spin 数量)
norms = [15279.4717, 15276.5560, 15277.7044, 15278.8215, 15279.2357, 
         15279.0916, 15277.4256, 15278.8801, 15278.3514, 15278.1166, 
         15278.3610, 15274.1688, 15277.0751, 15275.9344, 15275.6127, 
         15280.5963, 15281.7600, 15276.0051, 15274.7706, 15279.7130]

Mdots = [0.0]

# --- 2. 初始化 FITS 结构 ---
prihdr = fits.Header()
prihdr['INFO'] = "Lamppost corona intensity profile for an accretion disk of finite thickness"
prihdr['AUTHORS'] = 'BINGBING'
prihdr['VERSION'] = '1.2'
prihdu = fits.PrimaryHDU(header=prihdr)
aux = [prihdu]

# --- 3. 核心循环：处理物理参数 ---
k = 0  # 记录总的表格数量 (HDU index)

for q in Q_list:
    # 根据当前 Q 计算 spin 序列
    spin_max = np.sqrt(1 - q**2)
    spins = [0.0, 0.5*(spin_max-0.230), spin_max-0.230, spin_max-0.207, spin_max-0.185, 
             spin_max-0.165, spin_max-0.147, spin_max-0.129, spin_max-0.112, spin_max-0.097, 
             spin_max-0.081, spin_max-0.067, spin_max-0.053, spin_max-0.0265, spin_max-0.014, 
             spin_max-0.0018]

    # 将当前 Q 对应的 Spin 信息存入一个单独的 BinTable (可选)
    spin_col = fits.Column(name=f"spin_q_{q}", format="1E", array=np.array(spins))
    aux.append(fits.BinTableHDU.from_columns([spin_col]))

    for ii in range(len(spins)):
        spin = spins[ii]
        norm = norms[ii] if ii < len(norms) else norms[-1] # 防止索引溢出
        
        # 计算当前 spin 下的高度网格 (hgrid2)
        hmin0 = 1.099 * (1 + np.sqrt(1 - spin**2))
        hgrid2 = np.power(np.arange(100) / 99.0, 2.5) * (50.0 - hmin0) + hmin0

        Cols = []
        l = 0
        
        for h in hgrid2:
            # 构造文件名 (注意：%f 的精度需与磁盘文件名一致，建议使用 %.6f)
            filename = os.path.join(folder, "lp_%f_%f_%f.dat" % (spin, h,q))
            
            if os.path.exists(filename):
                try:
                    # 加载数据
                    rdisk, intensity, emis_delta, inc_delta = np.loadtxt(filename, unpack=True, skiprows=1)
                    
                    # 第一列存半径 r，后续列存物理量
                    if not Cols:
                        Cols.append(fits.Column(name='r_%d' % k, format='1E', array=rdisk))
                    
                    Cols.append(fits.Column(name='intensity_%d_%d' % (k, l), format='1E', array=intensity / norm))
                    Cols.append(fits.Column(name='emis_delta_%d_%d' % (k, l), format='1E', array=emis_delta))
                    Cols.append(fits.Column(name='inc_delta_%d_%d' % (k, l), format='1E', array=inc_delta))
                    
                    l += 1
                except Exception as e:
                    print(f"Error loading {filename}: {e}")
            else:
                # 如果找不到文件，打印提示（可选）
                print(f"Warning: File {filename} not found.")
                pass

        # 如果该 spin 下成功加载了数据，则创建一个新的 Table HDU
        if Cols:
            cols_defs = fits.ColDefs(Cols)
            tbhdu = fits.BinTableHDU.from_columns(cols_defs)
            aux.append(tbhdu)
            k += 1

# --- 4. 保存文件 ---
thdulist = fits.HDUList(aux)
try:
    thdulist.writeto(output_file, overwrite=True)
    print(f"Successfully saved to {output_file}")
except Exception as e:
    print(f"Failed to write FITS file: {e}")
finally:
    thdulist.close()

    import os
import numpy as np
from astropy.io import fits

folder = "data/"




prihdr = fits.Header()
prihdr['INFO'] = "Lamppost corona intensity profile for an accretion disk of finite thickness"
prihdr['AUTHORS'] = 'A. Abdikamalov, C. Bambi'
prihdr['COMMENTS'] = "Makes a file containing tabulated intensity values."
prihdu = fits.PrimaryHDU(header=prihdr)
aux=[prihdu]


files = os.listdir("data/")
spins = [-0.998, -0.75, -0.5, -0.25, 0.0, 0.2, 0.35, 0.5, 0.6, 0.69, 0.77, 0.8373257, 0.8689509, 0.8939505, 0.91543156, 0.9346257, 0.9521743, 0.9684631, 0.9837458, 0.9982]
rmins = [8.643391, 7.9241004, 7.1784782, 6.402528, 5.6693025, 5.067473, 4.520775, 4.034323, 3.6376512, 3.2537637, 2.8874276, 2.6138475, 2.442025, 2.2888975, 2.1422594, 1.9949745, 1.8398409, 1.6649805, 1.4352659, 1.2274946]
Mdots = [0, 0.05, 0.1, 0.2, 0.3]
norms = [15279.471701493161, 15276.556011507395, 15277.704447002663, 15278.821526714228, 15279.235704605322, 15279.091634291992, 15277.42562267709, 15278.880186568247, 15278.351395483163, 15278.116645095313, 15278.361014946237, 15274.168824563556, 15277.075106446524, 15275.934389911306, 15275.612738793216, 15280.596280508325, 15281.759979442828, 15276.005117296218, 15274.77055351251, 15279.713049336997]

cola = fits.Column(name="spin", format="1E", array = spins)
cola = fits.ColDefs([cola])
tbhdu = fits.BinTableHDU.from_columns(cola)
aux.append(tbhdu)

k = 0
# for spin in spins:
for ii in range(len(spins)):
	spin = spins[ii]
	norm = norms[ii]
	hmin1 = 1.7 * (1 + np.sqrt(1 - spin**2))
	hgrid1 = np.power(np.arange(40) / 39.0, 2) * (50.0 - hmin1) + hmin1
	hmin0 = 1.099 * (1 + np.sqrt(1 - spin**2))
	hgrid2 = np.power(np.arange(100) / 99.0,  2.5) * (50.0 - hmin0) + hmin0
	for mdot in Mdots:
	# for h in hgrid:
		Cols = []
		l = 0
		l2 = 0
		# for mdot in Mdots:
		for h in hgrid2:
			filename = "data/lp_%f_%f_%f.dat" % (spin, h, mdot)
			rdisk, intensity, emis_delta, inc_delta = np.loadtxt(filename, unpack=True, skiprows=1)
			if not len(Cols):
				Cols.append(fits.Column(name='r_%d' % k, format='1E',array=rdisk))
			Cols.append(fits.Column(name='intensity_%d_%d' % (k, l), format='1E',array=intensity / norm))
			Cols.append(fits.Column(name='emis_delta_%d_%d' % (k, l), format='1E',array=emis_delta))
			Cols.append(fits.Column(name='inc_delta_%d_%d' % (k, l), format='1E',array=inc_delta))
			l+=1


		cols = fits.ColDefs(Cols)
		tbhdu = fits.BinTableHDU.from_columns(cols)
		aux.append(tbhdu)
		k+=1


thdulist = fits.HDUList(aux)
thdulist.writeto('lp_thickness_v1.1b.fits')
thdulist.close()


