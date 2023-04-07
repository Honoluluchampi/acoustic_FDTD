import numpy as np
import matplotlib.pyplot as plt
import imageio.v2 as imageio
from IPython.display import display, Image
import os

class pressure_data:
  def __init__(self, xlen, ylen, data) :
    self.xlen = xlen
    self.ylen = ylen
    self.data = data

# ファイルを開いて、テキストデータを読み取ります
with open("output.txt") as f:
    # テキストデータを行ごとにリストに読み込みます
    data_lines = f.readlines()
    gif_name = data_lines[0]
    xsize = int(data_lines[1])
    ysize = int(data_lines[2])

# リストを numpy 二次元配列に変換します
data_array = np.zeros((len(data_lines), xsize, ysize))

for t in range(3, len(data_lines)):
    # コンマで区切った数値のリストを作ります
    line = data_lines[t]
    values = line.strip().split(",")
    # numpy 二次元配列に値を格納します
    for j in range(ysize):
      for i in range(xsize):
        data_array[t, i, j] = values[i + j * xsize]

data = pressure_data(xsize, ysize, data_array)

def save_pressure_gif(grid_data, name):
  # 各フレームの画像を生成し、GIFに追加
  frames = []
  dir_name = name
  anim_name = dir_name + "/animation.gif"
  if not os.path.exists(dir_name) : 
    os.mkdir(dir_name)

  print("creating gif file...")

  for i in range(0, len(grid_data.data)):
      grid = grid_data.data[i]
      
      # 生成した数値場を画像として保存
      filename = dir_name + '/frame' + str(i) + '.png'

      plt.imshow(grid, cmap='coolwarm', extent=(0, grid_data.xlen, 0, grid_data.ylen), 
                  vmin = -128, vmax = 128)
      plt.colorbar()
      plt.savefig(filename)
      plt.close()
      
      # 画像をGIFに追加
      frames.append(imageio.imread(filename).astype(np.uint8))

  # GIFを保存
  imageio.mimsave(anim_name, frames, duration=0.05)

  # 生成されたGIFファイルの読み込みと表示
  with open(anim_name,'rb') as f:
      display(Image(data=f.read(), format='gif'))


save_pressure_gif(data, gif_name)