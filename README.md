[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.18297314.svg)](https://doi.org/10.5281/zenodo.18297314)
# Chrono-Navier-Stokes-2D / Хроно-Навье-Стокса-2D

**Автор:** Максим Ковалевыч (maksymkovalevych)  
**Дата первой публикации:** январь 2026  
**Copyright © 2025–2026 Maksym Kovalevych. Все права защищены** (пока без лицензии — all rights reserved).

## О проекте / About the project

**Русский:**  
Chrono-NS — это гипотеза модификации уравнений Навье-Стокса, где собственное время (proper time τ) локально замедляется в зонах высокой завихрённости (vorticity ω):  

**Dτ/Dt = exp(-K₀ |ω|²)**  

Это создаёт эффект "хроновязкости" (chronoviscosity) — время становится адаптивной, саморегулирующейся средой, которая стабилизирует вихри и предотвращает blow-up (сингулярности).

Ключевые эффекты из 2D спектральной симуляции (Re ≈ 10⁴):  
- **Immersive Armor** («Иммерсивная броня») — вихри создают "временной барьер", защищая ядро.  
- **Statistical Safety Cut-off** — резкий обрыв тяжёлых хвостов в PDF завихрённости.  
- **Adaptive Energy Filtering** — спектр энергии гасит мелкие масштабы быстрее, усиливая обратный каскад.

Связь с реальностью: объясняет сверхдолговечные вихри, например Большое Красное Пятно Юпитера (350+ лет стабильности).

**English:**  
Chrono-NS is a hypothesis modifying Navier-Stokes equations where proper time τ dilates locally in high-vorticity zones:  

**Dτ/Dt = exp(-K₀ |ω|²)**  

This introduces "chronoviscosity" — time as an adaptive, self-regulating medium that stabilizes vortices and prevents blow-up.

Key effects from 2D spectral simulation (Re ≈ 10⁴):  
- **Immersive Armor** — vortex cores create a "temporal barrier" for protection.  
- **Statistical Safety Cut-off** — sharp exponential suppression of PDF tails.  
- **Adaptive Energy Filtering** — energy spectrum damps high-k noise, reinforces inverse cascade.

Real-world link: explains persistent coherent structures like Jupiter's Great Red Spot (350+ years).

## Файлы / Files
- `chrono_simulation.py` — основной скрипт симуляции (Standard vs Chrono-NS, дашборд).
- `Kovalevich_Logos_Chrono_Navier_Stokes_Simulation_v1.pdf` — подробный пост с визуалами и "THE LOGOS SOLUTION".

## Как запустить / How to run
1. Установи зависимости: `pip install numpy matplotlib scipy`
2. Запусти: `python chrono_simulation.py`
3. Получи дашборд сравнения Standard vs Chrono-NS.

## Лицензия / License
Пока **без лицензии** — все права защищены. Для использования/сотрудничества — обращайтесь.

Контакт: maksymkovalevych [at] gmail.com или через GitHub.

#ComputationalPhysics #NavierStokes #Turbulence #ChronoNS #KovalevichLogos
