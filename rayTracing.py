import numpy as np
from tqdm import tqdm

# Define a 3D vector class Vec3
# 3DベクトルクラスVec3を定義する
# 定义三维向量类 Vec3
class Vec3:
    # Initialize the 3D vector, default values are (0.0, 0.0, 0.0)
    # 3Dベクトルを初期化し、デフォルト値は(0.0, 0.0, 0.0)
    # 初始化三维向量，默认值为 (0.0, 0.0, 0.0)
    def __init__(self, e0=0.0, e1=0.0, e2=0.0):
        self.e = np.array([e0, e1, e2])

    # Get the x component of the vector
    # ベクトルのx成分を取得する
    # 获取向量的 x 分量
    def x(self):
        return self.e[0]

    # Get the y component of the vector
    # ベクトルのy成分を取得する
    # 获取向量的 y 分量
    def y(self):
        return self.e[1]

    # Get the z component of the vector
    # ベクトルのz成分を取得する
    # 获取向量的 z 分量
    def z(self):
        return self.e[2]

    # Define the negation operation of the vector
    # ベクトルの反転演算を定義する
    # 定义向量取反运算
    def __neg__(self):
        return Vec3(-self.e[0], -self.e[1], -self.e[2])

    # Define the element access by index
    # インデックスによる要素アクセスを定義する
    # 定义通过索引访问向量元素
    def __getitem__(self, i):
        return self.e[i]

    # Define the element setting by index
    # インデックスによる要素設定を定義する
    # 定义通过索引设置向量元素
    def __setitem__(self, i, value):
        self.e[i] = value

    # Define in-place addition of vectors
    # ベクトルのインプレース加算を定義する
    # 定义向量加法的原地运算
    def __iadd__(self, other):
        self.e += other.e
        return self

    # Define in-place multiplication of vector by scalar
    # スカラーによるベクトルのインプレース乗算を定義する
    # 定义向量乘标量的原地运算
    def __imul__(self, t):
        self.e *= t
        return self

    # Define in-place division of vector by scalar
    # スカラーによるベクトルのインプレース除算を定義する
    # 定义向量除以标量的原地运算
    def __itruediv__(self, t):
        self.e /= t
        return self

    # Define addition of vectors
    # ベクトルの加算を定義する
    # 定义向量加法
    def __add__(self, other):
        return Vec3(*(self.e + other.e))

    # Define subtraction of vectors
    # ベクトルの減算を定義する
    # 定义向量减法
    def __sub__(self, other):
        return Vec3(*(self.e - other.e))

    # Define multiplication of vector by scalar or vector
    # スカラーまたはベクトルによるベクトルの乗算を定義する
    # 定义向量与标量或向量相乘
    def __mul__(self, t):
        if isinstance(t, Vec3):
            return Vec3(*(self.e * t.e))
        return Vec3(*(self.e * t))

    # Define scalar multiplication by vector
    # スカラーによるベクトルの乗算を定義する
    # 定义标量与向量相乘
    def __rmul__(self, t):
        return self.__mul__(t)

    # Define division of vector by scalar
    # スカラーによるベクトルの除算を定義する
    # 定义向量除以标量
    def __truediv__(self, t):
        return Vec3(*(self.e / t))

    # Compute the length of the vector
    # ベクトルの長さを計算する
    # 计算向量长度
    def length(self):
        return np.linalg.norm(self.e)

    # Compute the squared length of the vector
    # ベクトルの長さの二乗を計算する
    # 计算向量长度的平方
    def length_squared(self):
        return np.dot(self.e, self.e)

# Normalize the vector
# ベクトルを正規化する
# 归一化向量
def unit_vector(v):
    return v / v.length()

# Define color as Vec3 type
# Vec3タイプとして色を定義する
# 定义颜色为 Vec3 类型
Color = Vec3
# Define point as Vec3 type
# Vec3タイプとして点を定義する
# 定义点为 Vec3 类型
Point3 = Vec3

# Define a ray class
# レイクラスを定義する
# 定义光线类
class Ray:
    # Initialize the ray, default origin and direction are (0.0, 0.0, 0.0)
    # レイを初期化し、デフォルトの原点と方向は(0.0, 0.0, 0.0)
    # 初始化光线，默认起点和方向均为 (0.0, 0.0, 0.0)
    def __init__(self, origin=Vec3(), direction=Vec3()):
        self.orig = origin
        self.dir = direction

    # Get the origin of the ray
    # レイの原点を取得する
    # 获取光线的起点
    def origin(self):
        return self.orig

    # Get the direction of the ray
    # レイの方向を取得する
    # 获取光线的方向
    def direction(self):
        return self.dir

    # Calculate the position of the ray at parameter t
    # パラメータtにおけるレイの位置を計算する
    # 计算光线在参数 t 处的位置
    def at(self, t):
        return self.orig + t * self.dir

# Define a function to check if a ray intersects a sphere
# レイが球と交差するかどうかを確認する関数を定義する
# 定义判断光线是否与球体相交的函数
def hit_sphere(center, radius, r):
    oc = r.origin() - center
    a = r.direction().length_squared()
    half_b = np.dot(oc.e, r.direction().e)
    c = oc.length_squared() - radius * radius
    discriminant = half_b * half_b - a * c
    if discriminant < 0:
        return -1
    else:
        return (- half_b - np.sqrt(discriminant)) / a

# Define a function to compute the color of a ray
# レイの色を計算する関数を定義する
# 定义计算光线颜色的函数
def ray_color(r):
    t = hit_sphere(Point3(0, 0, -1), 0.5, r)
    if t > 0.0:
        N = unit_vector(r.at(t) - Vec3(0, 0, -1))
        return 0.5 * Color(N.x() + 1, N.y() + 1, N.z() + 1)

    unit_direction = unit_vector(r.direction())
    a = 0.5 * (unit_direction.y() + 1.0)
    return (1.0 - a) * Color(1.0, 1.0, 1.0) + a * Color(0.5, 0.7, 1.0)

# Define a function to convert color values to a string
# 色の値を文字列に変換する関数を定義する
# 定义函数将颜色值转换为字符串
def write_color(pixel_color):
    r = pixel_color.x()
    g = pixel_color.y()
    b = pixel_color.z()

    # Convert the component values from [0,1] range to byte range [0,255]
    # [0,1]範囲の成分値をバイト範囲[0,255]に変換する
    # 将 [0,1] 范围的分量值转换为字节范围 [0,255]
    r_byte = int(255.999 * r)
    g_byte = int(255.999 * g)
    b_byte = int(255.999 * b)

    # Return the string representing the color components
    # 色成分を表す文字列を返す
    # 返回表示颜色分量的字符串
    return f"{r_byte} {g_byte} {b_byte}"

def main():
    # Image settings
    # 画像設定
    # 图像设置
    aspect_ratio = 16.0 / 9.0
    image_width = 1000
    image_height = int(image_width / aspect_ratio)
    image_height = max(image_height, 1)

    # Camera settings
    # カメラ設定
    # 相机设置
    focal_length = 1.0
    viewport_height = 2.0
    viewport_width = viewport_height * (image_width / image_height)
    camera_center = Point3(0, 0, 0)

    viewport_u = Vec3(viewport_width, 0, 0)
    viewport_v = Vec3(0, -viewport_height, 0)

    pixel_delta_u = viewport_u / image_width
    pixel_delta_v = viewport_v / image_height

    viewport_upper_left = camera_center - Vec3(0, 0, focal_length) - viewport_u / 2 - viewport_v / 2
    pixel00_loc = viewport_upper_left + 0.5 * (pixel_delta_u + pixel_delta_v)

    with open('rendered_image.ppm', 'w') as f:
        # Write the PPM header
        # PPMヘッダーを書く
        # 写入PPM头部
        f.write("P3\n")
        f.write(f"{image_width} {image_height}\n")
        f.write("255\n")

        # Render the image
        # 画像をレンダリングする
        # 渲染图像
        for j in tqdm(range(image_height), desc="Rendering"):
            for i in range(image_width):
                pixel_center = pixel00_loc + (i * pixel_delta_u) + (j * pixel_delta_v)
                ray_direction = pixel_center - camera_center
                r = Ray(camera_center, ray_direction)

                pixel_color = ray_color(r)
                f.write(write_color(pixel_color) + "\n")

if __name__ == "__main__":
    main()
