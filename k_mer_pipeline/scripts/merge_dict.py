import argparse

def process_sorted_files(file1_path, file2_path, output_path):
    with open(file1_path, 'r') as f1, open(file2_path, 'r') as f2, open(output_path, 'w') as out:
        print('start')
        # Итераторы с обработкой строк
        f1_iter = (line.strip().split() for line in f1)
        f2_iter = (line.strip().split() for line in f2)
        
        # Получаем первые значения из каждого файла
        parts1 = next(f1_iter, None)
        parts2 = next(f2_iter, None)

        print('cycle start')
        
        while parts1 is not None and parts2 is not None:
            if len(parts1) < 2 or len(parts2) < 2:  # Пропускаем некорректные строки
                parts2 = next(f2_iter, None)
                continue
            
            num1 = parts1[0]
            num2 = parts2[0]
            
            if int(num1) == int(num2):  # Совпадение - заменяем второе число
                
                if parts1[1] < parts2[1]:
                    out.write(f"{num2} {parts1[1]} {parts2[1]}\n")

                parts1 = next(f1_iter, None)
                parts2 = next(f2_iter, None)
                
            elif int(num1) < int(num2):  # Пропускаем меньшие значения из первого файла                
                parts1 = next(f1_iter, None)
            else:  # Записываем строку из второго файла без изменений
                out.write(f"{num2} 0 {parts2[1]}\n")
                parts2 = next(f2_iter, None)
        
        # Дозаписываем оставшиеся строки из второго файла
        while parts2 is not None:
            if len(parts2) >= 2:
                out.write(f"{parts2[0]} 0 {parts2[1]}\n")
            parts2 = next(f2_iter, None)

parser = argparse.ArgumentParser(description="Merge k-mers dictinaries")
parser.add_argument('-b', '--before', type=str, help = 'Dict of DSI (b - before)')
parser.add_argument('-a', '--after', type=str, help = 'Dict of DSII (b - after)')
parser.add_argument('-o', '--output', type=str, help = 'Output Merged Dict')

args = parser.parse_args()


# Пример использования
process_sorted_files(args.before, args.after, args.output) 
