#!/bin/bash

# Массив значений p (количество потоков)
THREAD_COUNTS=(1 2 3 4 5 6 7 8 9 10 20 40 60 80 100 120 140 160)

# Массив значений N (размеры входных данных)
SIZES=(10 100 1000 10000 100000)

# Количество повторений для каждого p и N
NUM_RUNS=5

# Исполняемый файл
PROGRAM="./var37_task_dynamic"

# Проверка существования исполняемого файла
if [ ! -x "$PROGRAM" ]; then
    echo "Исполняемый файл $PROGRAM не найден или не имеет прав на выполнение."
    exit 1
fi

# Файл для записи результатов
OUTPUT_FILE="results_var37_task_dynamic.txt"

# Инициализация файла результатов с заголовком
echo "N Threads Average_Time(s)" > "$OUTPUT_FILE"

# Цикл по каждому значению p (количество потоков)
for p in "${THREAD_COUNTS[@]}"
do
    echo "  Параллельные потоки: $p"

    # Инициализация ассоциативного массива для суммирования времени по N
    declare -A sum_time
    for N in "${SIZES[@]}"
    do
        sum_time[$N]=0.0
    done

    # Выполнение нескольких запусков для усреднения
    for run in $(seq 1 $NUM_RUNS)
    do
        # Установка количества потоков и запуск программы
        output=$(OMP_NUM_THREADS=$p "$PROGRAM")

        # Парсинг вывода программы для извлечения времени выполнения
        while IFS= read -r line
        do
            # Проверка соответствия строки формату "N=<значение> Time=<время>"
            if [[ $line =~ N=([0-9]+)\ Time=([0-9]+\.[0-9]+) ]]; then
                N="${BASH_REMATCH[1]}"
                time="${BASH_REMATCH[2]}"
                # Добавление времени к суммарному времени для текущего N
                sum_time[$N]=$(awk "BEGIN {printf \"%.6f\", ${sum_time[$N]} + $time}")
            fi
        done <<< "$output"
    done

    # Вычисление среднего времени для каждого N и запись в выходной файл
    for N in "${SIZES[@]}"
    do
        average_time=$(awk "BEGIN {printf \"%.6f\", ${sum_time[$N]} / $NUM_RUNS}")
        echo "$N $p $average_time" >> "$OUTPUT_FILE"
    done

    # Очистка ассоциативного массива для следующего p
    unset sum_time
done

echo "Результаты сохранены в файле: $OUTPUT_FILE"
