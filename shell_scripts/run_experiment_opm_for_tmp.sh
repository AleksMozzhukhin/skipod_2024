#!/bin/bash

# Массив значений p
THREAD_COUNTS=(1 2 3 4 5 6 7 8 9 10 20 40 60 80 100 120 140 160)

# Количество повторений для каждого p
NUM_RUNS=5

# Файл для записи результатов
OUTPUT_FILE="results_for_tmp.txt"

# Инициализация файла результатов
echo "Threads Average_Time(s)" > $OUTPUT_FILE

# Цикл по каждому значению p
for p in "${THREAD_COUNTS[@]}"
do
    echo "Running with OMP_NUM_THREADS=$p ..."
    total_time=0.0

    # Выполнение нескольких запусков
    for run in $(seq 1 $NUM_RUNS)
    do
        # Установка количества потоков и запуск программы
        output=$(OMP_NUM_THREADS=$p ./a.out)

        # Извлечение времени из вывода программы
        # Предполагается, что программа выводит строку вида "TOTAL_TIME: <время>"
        time=$(echo "$output" | grep "TOTAL_TIME" | awk '{print $2}')

        # Проверка корректности извлечения времени
        if [ -z "$time" ]; then
            echo "Failed to extract time for p=$p run=$run"
            continue
        fi

        # Суммирование времени
        total_time=$(awk "BEGIN {print $total_time + $time}")
    done

    # Вычисление среднего времени
    average_time=$(awk "BEGIN {print $total_time / $NUM_RUNS}")

    # Запись результатов в файл
    echo "$p $average_time" >> $OUTPUT_FILE
done

echo "Experiment completed. Results saved in $OUTPUT_FILE"
