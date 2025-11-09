'''
Программа моделирует движение нейтронов в 3Д области, состоящей из различных сред.

Каждый нейтрон рождается о случайной начальной позицией и направлением в пределах области источника. Область источника
выбирается из доступных областей при выполнении программы.
Для каждого нейтрона определяется длина свободного пробега (распределённая по экспоненциальному закону) 
исходя из характеристик(сечений) среды. Нейтрон движется прямолинейно до взаимодействия, после чего происходит одно из событий:
поглощение средой (нейтрон исчезает) или рассеяние под случайным углом относительно предыдущего направления.

Во время выполнения программы:
  *если необходимо, на экране отображается движение нейтронов в реальном времени;
  *параллельно в отдельном окне строится гистограмма отражающая поток на расстоянии r от центра цилиндров.
графикии строятся с помощью библиотеки MatPlotlib.
Расчёт потока нейтронов:
Поток рассчитывается как количество нейтронов, проходящих единицу площади за единицу времени.
(или, в данном случае, в единичном интервале радиуса в цилиндрических координатах). Для этого в каждой итерации моделирования
сохраняются текущие координаты всех нейтронов, после чего подсчитывается их число в каждом слое и делится на площадь слоя, строится гистограмма распределения.
Высота столбца гистограммы пропорциональна плотности нейтронов в данной области пространства,
что фактически отражает поток нейтронов на некотором расстоянии от оси (в конкретно данном примере, т.к. используется цилиндрическая геометрия).

Геометрия и входные данные:
Конкретно в данном примере рассматривается цилиндрическая геометрия, заданная через внешний OBJ-файл.
Этот файл содержит описание границ областей (среды и источника). OBJ файл может быть создан любым способом, но он не должен содержать
не замкнутых поверхностей(иначе программа не будет работать стабильно).

После каждого шага мёртвые (поглощённые) нейтроны заменяются новыми, это необходимо для поддержания
одинакового числа нейтронов, участвующих в симуляции(для стационарности процесса).
'''




# подключение библиотек и модулей                                                                                      
import trimesh                                            #3д объекты. импорт геометрии из файла
import numpy as np                                        # работа с массивами
import random as rnd                                      # генерация случайных величин
from math import pi, cos, sin, sqrt, log                  # необходимые математические функции 
import matplotlib.pyplot as plt                           # построение графика потока и отображение 3Д                 
from mpl_toolkits.mplot3d.art3d import Poly3DCollection   # работа с отображением 3Д вида                                        
from matplotlib.animation import FuncAnimation            # отвечает за анимацию                             


# функция спавна нейтронов в среде. В качестве входных параметров принимает среду, в которой необходимо создать нейтрон  и модуль его скорости
def spawn_neutron(environment, speed):

    # определяем координаты плоскостей, описывающих область
    minx, miny, minz = environment.mesh.bounds[0]
    maxx, maxy, maxz = environment.mesh.bounds[1]
    # случайным образом выбираем точку в описанном параллелепипеде до тех пор, пока выбранная точка не окажется в нужной области
    while True:
        # выбираем случайным образом координаты точки
        x = rnd.uniform(minx, maxx)
        y = rnd.uniform(miny, maxy)
        z = rnd.uniform(minz, maxz)
        # с помощью функции contains объекта библиотеки Triemesh проверяем, принадлежит ли точка области. если да, прекращаем цикл
        if environment.contains([x, y, z]):
            break

    # случайным образом выбираем направление движения
    phi = rnd.uniform(0, 2*pi)
    theta = rnd.uniform(0, pi)

    # рассчитываем компоненты скоростей
    v_x = speed * cos(phi) * sin(theta)
    v_y = speed * sin(phi) * sin(theta)
    v_z = speed * cos(theta)
    
    # функция возвращает объект класса Neutron  с рассчитаными параметрами
    return Neutron(x, y, z, v_x, v_y, v_z, environment)

#  класс описывающий среду. Содержит название среды, геометрическую область(mesh), сеченния поглощения и рассеяния, цвет и прозрачность(последние 2 параметра необходимы для отрисовки)
class Environment:
    def __init__(self, name, mesh, Sigma_a=0.0, Sigma_s=0.0, color='lightblue', opacity=0.5):
        self.name = name
        self.mesh = mesh
        # сечение поглощения
        self.Sigma_a = Sigma_a
        # сечение рассеяния
        self.Sigma_s = Sigma_s
        self.color = color
        self.opacity = opacity
    # вспомогательная функция. Определяет содержится ли точка point в данной среде
    def contains(self, point):
        # функция библиотеки Triemesh, проверяющая, содержит ли область точку
        return self.mesh.contains([point])[0]

# класс описывающий нейтрон. принимает начальные координаты, начальные компоненты скоростей и среду
class Neutron:
    def __init__(self, x, y, z, velocity_x, velocity_y, velocity_z, env):
        self.x = x
        self.y = y
        self.z = z
        self.velocity_x = velocity_x
        self.velocity_y = velocity_y
        self.velocity_z = velocity_z
        # модуль скорости нейтрона
        self.speed=sqrt(velocity_x**2+velocity_y**2+velocity_z**2)
        self.Enviroment = env
        # переменная отвечающая за путь пройденный нейтроном без взаимодействия
        self.free_way = 0
        # расчет длины свободного пробега в среде исходя из параметров среды
        self.lambda_free = -log(rnd.random()) / (self.Enviroment.Sigma_a + self.Enviroment.Sigma_s)
        # вспомогательные флаги, обозначающие, что
        self.alive = True       #нейтрон жив
        self.absorbed = False   #нейтрон поглощен

    # функция смещения нейтрона
    def update_Position(self):
        # временной шаг
        dt = 0.1
        # расчет новой позиции нейтрона
        new_x = self.x + self.velocity_x * dt
        new_y = self.y + self.velocity_y * dt
        new_z = self.z + self.velocity_z * dt

        # обновление пути нейтрона, пройденного без взаимодействия
        self.free_way += self.speed * dt

        # если путь пройденный без взаимодействия превосходит длину свободного пробега, разыгрываем взаимодействие и обнуляем его путь пройденный без взаимодействия
        if self.free_way >= self.lambda_free:
            self.interaction()
            self.free_way=0
        # проверка принадлежности нового положения нейтрона текущей среде
        if not self.Enviroment.contains([new_x, new_y, new_z]):
            # функция класса Neutron, возвращающая среду, в которой лежит новая точка
            self.determine_Enironment(new_x, new_y, new_z)
            # расчет длины нового свободного пробега (в новой среде) и обнуление пути, пройденного без взаимодействия
            self.lambda_free = -log(rnd.random()) / (self.Enviroment.Sigma_a + self.Enviroment.Sigma_s)
            self.free_way = 0
        
        # присвоение нейтрону новых координат
        self.x = new_x
        self.y = new_y
        self.z = new_z
    
    # функция класса Neutron. позволяет разыгрывать взаимодействия. Использует параметры самого объекта
    def interaction(self):
        # если Sigma_s / (Sigma_a + Sigma_s) > вероятность (в программе случайная величина от 0 до 1), то происходит рассивание нейтрона на случайный угол. Модуль скорости при этом не изменяется
        if self.Enviroment.Sigma_s / (self.Enviroment.Sigma_a + self.Enviroment.Sigma_s) > rnd.random():
            # задание случайного угла
            phi = rnd.uniform(0, 2*pi)
            theta = rnd.uniform(0, pi)
            # определение новых компонент скорости
            self.velocity_x = self.speed * cos(phi) * sin(theta)
            self.velocity_y = self.speed * sin(phi) * sin(theta)
            self.velocity_z = self.speed * cos(theta)
            # определение нового свободного пробега, обнуление уже пройденный путь
            self.lambda_free = -log(rnd.random()) / (self.Enviroment.Sigma_a + self.Enviroment.Sigma_s)
            self.free_way = 0
        else:
            # в другом случае нейтрон поглотился. Указываем, что нейтрон более не жив
            self.alive = False
            self.absorbed = True

    # вспомогательная функция, определяющая какой среде из списка(список ENVIRONMENTS_LIST-глобальная переменная) принадлежит нейтрон
    def determine_Enironment(self, new_x, new_y, new_z):
        # вспомогательный флаг, позволяющий определить, что новая среда определена корректно
        new=False
        for env in ENVIRONMENTS_LIST:
            # проверяем каждую среду на наличие точки
            if env.contains([new_x, new_y, new_z]):
                # если нашли, присваиваем нейтрону новую среду и меняем специальный флаг
                self.Enviroment = env
                new=True
        # если новая среда не найдена, значит нейтрон улетел за пределы области симуляции. Ставим отметку о том, что нейтрон мертв
        if not new: self.alive = False

#Загружаем объекты из файла .OBJ, в дополнительно необходим файл .mtl, содержащий данные о материалах
# в данном примере используется тонкий цилиндрический источник в большой среде.
scene = trimesh.load("source_line.obj", force='scene')

# глобальная переменная. Массив всех сред
ENVIRONMENTS_LIST = []

# вспомогательная перменная. Необходима для масштабирования 3д моделей из obj файла
scale_factor = 1000  

# из scene импортируются названия материалов и областей
for name, geom in scene.geometry.items():
    # масштабируем объекты
    geom.apply_scale(scale_factor)
    color = 'lightblue'

    # поиск известных сред по названиям и присвоение им нужных параметров(сечений, цвета)
    if 'environment' in name.lower():
        color = 'blue'
        Sigma_a = 0.003
        Sigma_s = 0.1
    elif 'source' in name.lower():
        color = 'green'
        Sigma_a = 0.000000001
        Sigma_s = 0.000000001
    elif 'reflector' in name.lower():
        color = 'gray'
        Sigma_a = 0.001
        Sigma_s = 1.5
    # добавление среды в список сред
    ENVIRONMENTS_LIST.append(Environment(name, geom, Sigma_a=Sigma_a, Sigma_s=Sigma_s, color=color, opacity=0.3))


# ///////////////////////////////////////////////////////////////////////////////////////////////

# количество нейтронов
n_neutrons = 1000
enable_3d_view = True  #  поставить False, если Вы не хотите увидеть очень красивую 3D-анимацию

# массив для хранения нейтронов
neutrons = []

# ///////////////////////////////////////////////////////////////////////////////////////////////

# блок определения источника нейтронов. сначала источником делается случайная среда(сделано во избежания неопределенностей)
source_environment = rnd.choice(ENVIRONMENTS_LIST)
# поиск в списке среды с названием "source". Найденная среда является источником
for envr in ENVIRONMENTS_LIST:
    if envr.name == 'source':
        source_environment = envr

# создание n_neutrons нейтронов в источнике
for _ in range(n_neutrons):
    neutrons.append(spawn_neutron(source_environment, speed=5))

# сохранение всех плоскостей в массив. необходимо для дальнейшей отрисовки
all_vertices = np.vstack([env.mesh.vertices for env in ENVIRONMENTS_LIST])


# вспомогательная переменная. необходима для расчета потока в примере с цилиндрическим тонким источником
H = all_vertices[:, 2].max() - all_vertices[:, 2].min()

#блок визуализации
# если необходимо отобразить 3Д вид
if enable_3d_view:
    # создается объект fig из MatPlotlib. это окно в котором отобразится 3д вид и гистограмма потока
    fig = plt.figure(figsize=(13, 5))
    ax3d = fig.add_subplot(121, projection='3d')
    ax_hist = fig.add_subplot(122)
else:
    # создается объект fig из MatPlotlib. это окно в котором отобразится гистограмма потока
    fig = plt.figure(figsize=(6, 5))
    ax3d = None
    ax_hist = fig.add_subplot(111)

# если включен 3Д вид, происходит преобразование их объектов библиотеки Triemesh в объекты библиотеки MatPlotLib
if enable_3d_view:
    for env in ENVIRONMENTS_LIST:
        faces = env.mesh.faces
        vertices = env.mesh.vertices
        mesh = [[vertices[i] for i in face] for face in faces]
        poly = Poly3DCollection(mesh, alpha=env.opacity, facecolor=env.color, edgecolor='k', linewidths=0.1)
        ax3d.add_collection3d(poly)

    # создание массива с координатами всех нейтронов
    neutron_scatter = ax3d.scatter(
        [n.x for n in neutrons],
        [n.y for n in neutrons],
        [n.z for n in neutrons],
        color='blue', s=5
    )

    # формирование 3д вида их всех плоскостей 
    ax3d.set_xlim(all_vertices[:, 0].min(), all_vertices[:, 0].max())
    ax3d.set_ylim(all_vertices[:, 1].min(), all_vertices[:, 1].max())
    ax3d.set_zlim(all_vertices[:, 2].min(), all_vertices[:, 2].max())
    ax3d.set_xlabel('X')
    ax3d.set_ylabel('Y')
    ax3d.set_zlabel('Z')
    ax3d.set_title("3Д вид")

#формирование гистограммы
# задаются столбцы гистограммы: np.linspace(от начальной координаты, до конечной, на сколько частей делить)
bins = np.linspace(0, 20, 200)
hist_data = ax_hist.bar(bins[:-1], np.zeros(len(bins)-1), width=np.diff(bins), align='edge', color='purple', alpha=0.7)
ax_hist.set_xlim(0, bins[-1])
ax_hist.set_ylim(0, n_neutrons)
ax_hist.set_xlabel("Расстояние от центра")
ax_hist.set_ylabel("Плотность нейтронов")
ax_hist.set_title("Нейтронный поток")


# функция отвечающая за 1 итерацию симуляции. Требует 1 аргумента на вход. Аргумент не используется в самой функции
def update(frame):
    # вспомогательный счетчик мертвых нетронов
    dead=0
    # для каждого живого нейтрона обновляется его позиция. если нейтрон поглощен или вылетел за пределы области симуляции, счетчик "мертвых" нейтронов увеличивается на 1
    for p in neutrons:
        if p.alive: p.update_Position()
        if not p.alive: dead += 1 

    # удаление из массива всех "мертвых" нейтронов
    neutrons[:] = [n for n in neutrons if n.alive]

    # восполнение потерей нейтронов. В источнике создаются нейтроны в количестве равном количеству потеряных нейтронов
    for _ in range(dead):
        neutrons.append(spawn_neutron(source_environment, speed=5))

    # создание массива с координатами всех живых нейтронов вне источника. Необходимо для расчета потока
    pts = np.array([[n.x, n.y, n.z] for n in neutrons if n.alive and n.Enviroment != source_environment])

    # блок расчета потока нейтронов. Для примера с цилиндрической геометрией
    '''
    Порядок расчета потока:
    Берется радиус, на котором необходимо расчитать поток. Делаются отступы к и от центра. Считается число нейтронов
    в получившемся слое. Количество нетйронов в слое делится на площадь цилиндрического слоя начального радиуса.

    ПРИМЕЧАНИЕ: в данной программе данная схема реализована следующим образом:
    Из массива bins, содержащего начальную и конечную координаты столбцов формируются массивы всех левых 
    и всех правых границ. создается массив всех средних радиусов.
    '''
    # если число частиц не нулевое
    if len(pts) > 0:

        #Перевод координат (x y) в радиус
        r_cyl = np.sqrt(pts[:, 0]**2 + pts[:, 1]**2)

        # расчет радиусов левой и правой границ
        # массив всех ПРАВЫХ границ
        r_left = bins[:-1]
        # массив всех ЛЕВЫХ границ
        r_right = bins[1:]
        # массив всех СРЕДНИХ радиусов
        r_mid = 0.5 * (r_left + r_right)

        # расчет площади цилиндрической оболочки
        S = 2 * np.pi * r_mid * H
        # во избежание деления на ноль
        S[S == 0] = 1e-10

        # функция np.histogram(данные, столбцы) возвращает массив с количеством элементов в каждом бине (counts) и массив границ (_). Массив границ далее не используется
        counts, _ = np.histogram(r_cyl, bins)

        # расчет потока как отношение числа нейтронов на площадь умноженное на масштабирующий множитель
        flux_density = counts / S *1e3 # нейтроны на м2

        # функция zip объединяет массивы hist_data и flux_density в один двумерный (первый из hist_data с первым из flux_density; воторой из hist_data со вторым из flux_density и т.д.)
        # в цикле задается значение переменной rect из массива hist_data, h из flux_density. Таким образом каждому прямоугольнику соответсвует его новая высота
        for rect, h in zip(hist_data, flux_density):
            # установка для каждого прямоугольника его высоты
            rect.set_height(h)

    # если число частиц нулевое, у всех прямоугольников высота нулевая
    else:
        for rect in hist_data:
            rect.set_height(0)

    # установка диапазона по оси ординат для гистограммы
    ax_hist.set_ylim(0, 10)

    # если включено отображение 3Д вида и число частиц не 0
    if enable_3d_view and len(pts) > 0:
        # обновление всех нейтронов на 3Д виде
        neutron_scatter._offsets3d = (pts[:, 0], pts[:, 1], pts[:, 2])

    # print(len(neutrons))

    # возвращение сворфмированных массивов. если 3Д включено, и данные для гистограммы и для 3Д. Иначе только для гистограммы
    return (*hist_data,) if not enable_3d_view else (*hist_data, neutron_scatter)


# функция MatPlotLib, вызывающая функцию update. Именно она отвечает за проведение симуляции и отображения 3Д вида и гистограммы
anim = FuncAnimation(fig, update, frames=range(300), interval=100, blit=False)
# отображение графиков.
plt.show()

