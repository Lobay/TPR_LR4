from copy import deepcopy


INF = 1e100  # бесконечность
EPSILON = 1e-10  # точность
MIN = 0  # режим минимизации
MAX = 1  # режим максимизации


# получение вещественной части числа
def GetRealPart(x: float) -> float:
    if x > 0:
        return x - int(x)

    return int(-x) + 1 + x


class SimplexSolve:
    def __init__(self, x, f):
        self.x = x
        self.f = f

    # вывод решения
    def Print(self):
        for i, xi in enumerate(self.x):
            print("x" + str(i + 1), "=", str(xi) + ", ", end='')

        print("F:", self.f)

    # получение индекса вещественного решения
    def GetRealIndex(self):
        for i, xi in enumerate(self.x):
            if abs(xi - int(xi)) > EPSILON:  # если не совпадает с целой частью
                return i  # значит нашли

        return -1  # нет такого

    # получение индекса вещественного решения с максимальной вещественной частью
    def GetRealMaxIndex(self):
        imax = -1

        for i, xi in enumerate(self.x):
            if xi != int(xi) and (imax == -1 or GetRealPart(self.x[imax]) < GetRealPart(xi)):
                imax = i

        return imax

    # получение значения решения
    def GetX(self, index):
        return self.x[index]

    # получение значения функции
    def GetF(self):
        return self.f


class SimplexTable:
    # конструктор симплекс-таблицы
    def __init__(self, a, b, c, mode):
        self.a = deepcopy(a)
        self.b = deepcopy(b)
        self.c = deepcopy(c)
        self.mode = mode
        self.basis = []

        self.m = len(a)  # запоминаем количество ограничений
        self.n = len(a[0])  # запоминаем количество переменных

        for i in range(self.m):
            for j in range(self.m):
                self.a[i].append(i == j)  # добавляем единичную матрицу
            
            self.c.append(0)  # добавляем нули в вектор функции
            self.basis.append(i + self.n)

        self.c.append(0)
        self.deltas = [0 for _ in range(self.n + self.m + 1)]  # дельты
        self.UpdateDeltas()

    # получение количества основных переменных
    def GetMainVariablesCount(self):
        return self.n

    # получение количества переменных
    def GetVariablesCount(self):
        return self.n + self.m

    # получение количества ограничений
    def GetRestrictionCount(self):
        return self.m

    def Print(self):
        print("| basis |", end='')
        for i in range(self.n + self.m):
            print("        x" + str(i + 1), "|", end='')
        print("         b |")

        for i in range(self.m):
            print("|    x" + str(self.basis[i] + 1), "|", end='')

            for j in range(self.n + self.m):
                print(" %9.5f |" % self.a[i][j], end='')

            print(" %9.5f |" % self.b[i])

        print("|     F |", end='')

        for i in range(self.n + self.m):
            print(" %9.5f |" % (self.deltas[i] if self.mode == MIN else -self.deltas[i]), end='')

        print(" %9.5f |" % (self.deltas[self.n + self.m] if self.mode == MIN else -self.deltas[self.n + self.m]))

    # вывод задачи
    def PrintTask(self):
        print("F = ", end='')

        for i in range(self.n):
            if self.c[i] >= 0 and i > 0:
                print("+", end='')

            print(str(self.c[i]) + "x" + str(i + 1), end='')

        print("")

        for i in range(self.m):
            for j in range(self.n):
                if self.a[i][j] >= 0 and j > 0:
                    print("+", end='')

                print(str(self.a[i][j]) + "x" + str(j + 1), end=' ')

            print("<=", self.b[i])

    # деление строки на число
    def DivideRow(self, row, value):
        for i in range(self.n + self.m):
            self.a[row][i] /= value

        self.b[row] /= value

    # вычитание строки из другой строки, умноженной на число
    def SubstractRow(self, row1, row2, value):
        for i in range(self.n + self.m):
            self.a[row1][i] -= self.a[row2][i] * value

        self.b[row1] -= self.b[row2] * value

    # исключение Гаусса
    def Gauss(self, row, column):
        self.DivideRow(row, self.a[row][column])  # делим строку на опорный элемент

        for i in range(self.m):
            if i != row:
                self.SubstractRow(i, row, self.a[i][column])

        self.basis[row] = column  # делаем выбранный элемент базисным

    # расчёт дельт
    def UpdateDeltas(self):
        for i in range(self.n + self.m + 1):
            self.deltas[i] = -self.c[i]

            for j in range(self.m):
                self.deltas[i] += self.c[self.basis[j]] * (self.b[j] if i == self.n + self.m else self.a[j][i])

    # получение дельты
    def GetDelta(self, index):
        return self.deltas[index]

    # получение значения матрицы
    def GetValue(self, row, column):
        return self.a[row][column]

    # получение свободного коэффициента
    def GetB(self, index):
        return self.b[index]

    # получение базисной переменной
    def GetBasis(self, index) -> int:
        return self.basis[index]

    # получение симплекс отношения
    def GetRelation(self, row, column):
        if abs(self.a[row][column]) < EPSILON:
            return INF
     
        if self.b[row] >= 0 and self.a[row][column] < 0:
            return INF
            
        return self.b[row] / self.a[row][column]

    # получение значения функции
    def GetF(self, x):
        f = 0

        for i in range(self.n):
            f += self.c[i] * x[i]  # считаем значение функции

        return f

    # получение решения
    def GetSolve(self):
        x = [0 for _ in range(self.n)]

        # заполняем из базисных переменных
        for i in range(self.m):
            if self.basis[i] < self.n:
                x[self.basis[i]] = self.b[i]

        return SimplexSolve(x, self.deltas[self.n + self.m])  # возвращаем решение

    # проверка ограничения
    def CheckRestriction(self, x, index):
        s = 0

        for i in range(self.n):
            s += self.a[index][i] * x[i]  # считаем линейную комбинацию

        return s <= self.b[index]  # проверяем ограничение

    def IsBasis(self, index):
        return index in self.basis

    def AddBasis(self, index):
        self.basis.append(index)

    def AddRestriction(self, b):
        for row in self.a:
            row.append(0)

        self.a.append([0 for _ in range(self.n + self.m + 1)])
        self.b.append(b)
        self.c.append(0)
        self.deltas.append(0)
        self.m += 1
