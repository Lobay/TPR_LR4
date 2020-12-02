from copy import deepcopy
from SimplexTable import SimplexSolve, SimplexTable, INF, EPSILON, MIN, MAX, GetRealPart


class Simplex:
    def __init__(self, a, b, c, mode):
        self.table = SimplexTable(a, b, c, mode)
        self.mode = mode

        self.taskA = deepcopy(a)
        self.taskB = deepcopy(b)
        self.taskC = deepcopy(c)

    # получение строки с максимальным по модулю отрицательным значением b
    def get_row_with_negative_b(self):
        index = -1

        for i in range(self.table.GetRestrictionCount()):
            if (self.table.GetB(i) < 0) and (index == -1 or self.table.GetB(i) < self.table.GetB(index)):
                index = i

        return index

    # получение столбца с максимальным по модулю элементом в строке
    def get_column_with_negative_b(self, row):
        index = -1

        for i in range(self.table.GetVariablesCount()):
            if self.table.GetValue(row, i) >= 0:
                continue

            if index == -1 or abs(self.table.GetValue(row, i)) > abs(self.table.GetValue(row, index)):
                index = i

        return index

    # удаление отрицательных элементов в b
    def exclude_negative_b(self):
        while True:
            row = self.get_row_with_negative_b()  # получаем строку с отрицательным b

            if row == -1:
                return True

            column = self.get_column_with_negative_b(row)  # получение столбца

            print("Remove negative b at row", row + 1)
            self.table.Print()  # выводим таблицу

            if column == -1: # если не нашли
                return False  # значит нельзя избавиться

            self.table.Gauss(row, column)  # выполняем исключение Гауса

    # проверка плана на оптимальность
    def is_optimal(self):
        for i in range(self.table.GetVariablesCount()):
            if self.mode == MAX and self.table.GetDelta(i) < 0:
                return False

            if self.mode == MIN and self.table.GetDelta(i) > 0:
                return False

        return True  # план оптимален

    # подходит ли решение по условию
    def is_solve(self, x):
        for i in range(self.table.GetRestrictionCount()):
            if not self.table.CheckRestriction(x, i): # если не выполнено ограничение
                return False  # то решение не подходит

        return True  # решение подходит

    # расчёт симплекс отношений
    def get_relations(self, column):
        q = []

        for row in range(self.table.GetRestrictionCount()):
            q.append(self.table.GetRelation(row, column))

        return q

    # получение разрешающего столбца
    def get_column(self):
        column = 0

        for i in range(self.table.GetVariablesCount()):
            if self.mode == MAX and self.table.GetDelta(i) <= self.table.GetDelta(column):
                column = i
            elif self.mode == MIN and self.table.GetDelta(i) >= self.table.GetDelta(column):
                column = i

        return column

    # получение разрешающей строки
    def get_row(self, q):
        row = -1

        for i in range(self.table.GetRestrictionCount()):
            if q[i] != INF and (row == -1 or q[i] <= q[row]):
                row = i

        return row

    # решение задачи
    def solve(self):
        if not self.exclude_negative_b():
            print("solve does not exists\n")
            return False

        iteration = 1

        while True:
            self.table.UpdateDeltas()  # расчитываем дельты

            print("\nIteration", iteration)
            iteration += 1
            self.table.Print()

            # проверка плана на оптимальность, если оптимален, то решение найдено
            if self.is_optimal():
                return True  # решение есть

            column = self.get_column()  # получаем разрешающий столбец
            q = self.get_relations(column)  # рассчитываем симлекс-отношения
            row = self.get_row(q)  # получаем разрешающую строку

            # если нет разрешающей строки, то решения нет
            if row == -1:
                print("solve does not exist\n")
                return False

            self.table.Gauss(row, column)  # выполняем исключение Гауса

    def solve_gomory(self):
        print("Solve task")
        self.table.PrintTask()

        if not self.solve():  # если решение вещественной задачи не было найдено
            return None  # то и решения нет

        solve = self.table.GetSolve()  # получаем решение
        solve.Print()
        print("")

        real_index = solve.GetRealMaxIndex()  # ищем вещественное число с максимальной дробной частью

        if real_index == -1:  # если решение целочисленное
            return solve  # возвраащем решение

        # ищем строку с этой базисной переменной
        index = 0
        while self.table.GetBasis(index) != real_index:
            index += 1

        self.table.AddRestriction(-GetRealPart(solve.GetX(real_index)))

        for i in range(self.table.GetVariablesCount()):
            self.table.a[self.table.GetRestrictionCount() - 1][i] = 0 if self.table.IsBasis(i) else -GetRealPart(self.table.a[index][i])

        self.table.a[-1][self.table.GetVariablesCount() - 1] = 1
        self.table.AddBasis(self.table.GetVariablesCount() - 1)

        print("Add GOMORY restriction")
        self.solve_gomory()

    # получение целочисленных решений
    def solve_integer_branches_and_borders(self, depth=0):
        print("Start solving task:")
        self.table.PrintTask()

        # если решение не было найдено, то добавляем пустое решение
        if not self.solve():
            return []

        solve = self.table.GetSolve()  # получаем решение
        solve.Print()
        print("")

        real_index = solve.GetRealIndex()  # ищем вещественное решение

        # если решение содержало только целые числа, то возвращаем решение
        if real_index == -1:
            return [solve]

        b = int(solve.GetX(real_index))  # получаем новое условие

        # готовим новые условия
        a1 = deepcopy(self.taskA)
        a2 = deepcopy(self.taskA)

        b1 = deepcopy(self.taskB)
        b2 = deepcopy(self.taskB)

        n = self.table.GetMainVariablesCount()
        m = self.table.GetRestrictionCount()

        # разбиваем задачу на 2 ветки решения
        a1.append([0 for _ in range(n)])
        a1[m][real_index] = 1
        b1.append(b)
        simplex1 = Simplex(a1, b1, deepcopy(self.taskC), self.mode)

        a2.append([0 for _ in range(n)])
        a2[m][real_index] = -1
        b2.append(-b - 1)
        simplex2 = Simplex(a2, b2, deepcopy(self.taskC), self.mode)

        print("Divide for x" + str(real_index + 1), "<=", b, "and x" + str(real_index + 1), ">=", b + 1)

        # запускаем решение для каждой ветки
        solve1 = simplex1.solve_integer_branches_and_borders(depth + 1)
        solve2 = simplex2.solve_integer_branches_and_borders(depth + 1)

        return solve1 + solve2  # возвращаем решения

    # поиск решений методом грубой силы
    def solve_integer_bruteforce(self, n_max):
        solves = []
        print("\nBruteforce solves: ");

        for x1 in range(n_max):
            for x2 in range(n_max):
                for x3 in range(n_max):
                    x = [x1, x2, x3 ]

                    if not self.is_solve(x):  # если решение не подходит под ограничения
                        continue  # то пропускаем

                    solve = SimplexSolve(x, self.table.GetF(x))
                    solve.Print()
                    solves.append(solve)  # добавляем решение

        return solves  # возвраащем решения

    # поиск лучшего из решений
    def find_best_solve(self, solves):
        index = 0

        for i, solve in enumerate(solves):
            if self.mode == MAX and solves[i].GetF() > solves[index].GetF():
                index = i
            elif self.mode == MIN and solves[i].GetF() < solves[index].GetF():
                index = i

        print("Best solve: ")
        solves[index].Print()
