import pandas as pd
import matplotlib.pyplot as plt

# Чтение данных из файла с указанием кодировки
data = pd.read_csv('result.txt', sep='\t', encoding='windows-1251')  # Попробуйте 'ISO-8859-1' если это не сработает

# Преобразование типов данных
data['Параметр_N'] = data['Параметр_N'].astype(int)
data['Количество_нитей'] = data['Количество_нитей'].astype(int)
data['Время_выполнения'] = data['Время_выполнения'].astype(float)

# Создание графика
plt.figure(figsize=(12, 6))

# Построение графиков для каждого количества нитей
for threads in data['Количество_нитей'].unique():
    subset = data[data['Количество_нитей'] == threads]
    plt.plot(subset['Параметр_N'], subset['Время_выполнения'], marker='o', label=f'Нитей: {threads}')

# Настройка графика
plt.title('Зависимость времени выполнения от параметра N и количества нитей')
plt.xlabel('Параметр N')
plt.ylabel('Время выполнения (секунды)')
plt.xscale('log')  # Логарифмическая шкала по оси X
plt.yscale('log')  # Логарифмическая шкала по оси Y
plt.legend()
plt.grid(True)
plt.tight_layout()

# Показать график
plt.show()