# 📘 Полное руководство пользователя и техническая документация

### Программа: Генерация и анализ двумерных гауссиан

---

## 🧭 Содержание

1. Введение
2. Назначение программы
3. Архитектура программы
4. Компиляция и запуск
5. Конфигурационный файл `config.txt`
6. Команды и синтаксис batch-файла
7. Подробное описание каждой команды
8. Алгоритмы: генерация, кластеризация, "волна"
9. Работа с изображениями BMP
10. Взаимодействие с GnuPlot
11. Примеры batch-файлов
12. Сценарии использования (use-cases)
13. Возможные ошибки и способы их устранения
14. Логирование и диагностика
15. Расширение программы
16. Лицензия и авторство

---

## 1. 📌 Введение

Программа написана на языке C и предназначена для генерации двумерных гауссовых распределений с последующей визуализацией и анализом. Она позволяет строить поле, накладывать на него гауссы, выполнять пороговые срезы, выделять компоненты, проводить кластеризацию и визуализировать результат в форматах BMP и GnuPlot.

---

## 2. 🎯 Назначение программы

Программа решает следующие задачи:

- Симуляция распределения плотности (например, интенсивности сигналов, тепловых карт и др.)
- Предобработка изображений для анализа компонент
- Простая кластеризация на основе KMeans
- Выделение связных компонент на основе алгоритма «волна»

---

## 3. 🧱 Архитектура программы

Файл `test.c` представляет собой консольное приложение. Ключевые модули:

- **main()** — основной цикл обработки
- **инициализация параметров** — через config-файл
- **парсинг batch-файлов**
- **реализация команд: INI, G, GEN, BIN, BMP, WAVE, KMEANS, RBMP**
- **визуализация** через текстовые файлы или BMP

---

## 4. ⚙️ Компиляция и запуск

```bash
gcc test.c -o gaussgen.exe
./gaussgen.exe
```

---

## 5. 🧾 Конфигурационный файл `config.txt`

Пример содержимого:

```txt
Log_Server = ON //Логирование на сервере
Log_Server_filename = log_server.txt // Файл логирования сервера
Log_Client = ON // Логирование на клиенте
Log_Client_filename = log_client.txt // Файл логирования клиента
Batch_pause = OFF // Пауза после каждой команды для отладки
Batch_file = ON // Командный файл
Batch_filename = Batch120.txt // Название командного файла
Field_DEF = 120 100 // Стандартные размеры поля
Gauss_DEF = 150 50 50 10 10 // Стандартный гауз для добавления
Noize = 11 // Компоненты до какого размера не пропускать
```

---

## 6. 🧑‍💻 Команды batch-файла

Каждая строка — отдельная команда. Пример:
```txt
INI 120 100
G 150 50 50 10 10
GEN
BIN 10
GNU
BMP result.bmp
WAVE
KMEANS 3 20
```

---

## 7. 🔍 Подробное описание команд

| Команда            | Описание |
|--------------------|----------|
| `INI X Y`          | Инициализация поля X на Y |
| `G A B C D E`      | Добавить гаусс с амплитудой A, центром B, C и сигмами D, E |
| `GEN`              | Генерация поля |
| `BIN H`            | Бинаризация по высоте H |
| `GNU`              | Создание файла для GnuPlot |
| `BMP fname.bmp`    | Генерация BMP |
| `WAVE`             | Запуск "волны" |
| `KMEANS K P`       | Кластеризация с K центрами и P точками |
| `RBMP file.bmp`    | Загрузка BMP-файла |

---

## 8. 🧠 Алгоритмы

### Генерация гаусса
```
f(x, y) = A * exp(-((x - x0)^2 / (2 * sigma_x^2) + (y - y0)^2 / (2 * sigma_y^2)))
```

### "Волна" — поиск связных компонент (аналог flood-fill)

### KMeans — итеративная кластеризация с двумя уровнями центроид

---

## 9. 🖼 Работа с BMP

- Генерация BMP по текущему состоянию поля
- Визуализация компонент, кластеров, треугольников и путей

---

## 10. 📈 Взаимодействие с GnuPlot

Программа создаёт файлы `.txt` для 3D-визуализации:
```gnuplot
splot "gnuplot_plot.txt" with lines
```

---

## 11. 🧪 Примеры batch-файлов

Сценарий генерации поля и анализа:
```txt
INI 200 150
G 100 50 50 10 10
G 80 130 70 15 15
GEN
BIN 10
WAVE
BMP output.bmp
KMEANS 3 10
GNU
```

---

## 12. 📂 Сценарии использования

- Моделирование сигналов
- Обработка тепловых карт
- Тестирование кластеризации и компонентного анализа
- Обучающие цели

---

## 13. 🧯 Возможные ошибки

| Ошибка                       | Решение                           |
|-----------------------------|-----------------------------------|
| "Field not initialized"     | Добавьте команду `INI` в batch    |
| "No gausses to apply"       | Убедитесь, что вы вызвали `G`     |
| "Invalid BMP file"          | Проверьте формат изображения      |

---

## 14. 🗃 Логирование

Включается через `config.txt`. Генерирует:

- `log_server.txt` — служебные логи
- `log_client.txt` — действия пользователя

---

## 15. 🧩 Расширение программы

Можно добавить:
- Алгоритмы машинного обучения
- Улучшение фильтрации изображений
- Поддержку других форматов ввода/вывода

---

## 16. 👨‍💻 Авторство

**Коротков Павел**, 212 группа.
