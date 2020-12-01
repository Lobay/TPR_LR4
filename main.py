from SimplexSolve import Simplex, MAX, MIN


def main():
    c = [1, 3, 8]
    a = [
        [1, 1, 1],
        [1, 1, 0],
        [0, 0.5, 2]
    ]

    b = [7, 2, 4]

    # # метод ветвей и границ
    # branches = Simplex(a, b, c, MAX)
    # branchesSolves = branches.solve_integer_branches_and_borders()  # ищем решения методом ветвей и границ
    # branches.find_best_solve(branchesSolves)  # находим лучшее решение
    #
    # # метод перебора
    # xmax = 20  # максимальное значения для перебора
    # bruteforce = Simplex(a, b, c, MAX)
    # bruteforceSolves = bruteforce.solve_integer_bruteforce(xmax)  # ищем решения методом перебора
    # bruteforce.find_best_solve(bruteforceSolves)  # находим лучшее решение

    gomory = Simplex(a, b, c, MAX)
    gomory.solve_gomory()


if __name__ == '__main__':
    main()
