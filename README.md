# sph-toy-star

![SPH Star GIF](sph.GIF "N = 200")

**[ссылка](https://vk.com/away.php?to=https%3A%2F%2Fwww.youtube.com%2Fwatch%3Fv%3Dbs45NErLxuA&cc_key=)** на видео с каким-никаким описанием

## как собрать проект?

1. сборка происходит с помощью утилиты cmake, так что cmake должен быть установлен :>
2. в терминале переходим в папку, в которой хотим оставить проект.
3. пишем: `git clone git@github.com:Shulliikkk/sph-toy-star.git`
4. `cd sph-toy-star`
5. `cmake -S . -B build` -- эта команда запустит cmake и скажет ему, где расположен сам проект и в какой папке его собирать 
6. `cmake --build build` -- собираем проект
7. `build/SPH-STAR` -- запуск программы ("SPH-STAR" -- всего лишь исполняемый файл) 

