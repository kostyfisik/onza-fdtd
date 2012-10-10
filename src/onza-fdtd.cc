// Using Doxygen 1.8.0 (with Markdown)
///
/// @file   onza-fdtd.cc
/// @author Ladutenko Konstantin <kostyfisik at gmail (.) com>
/// @copyright 2012 Ladutenko Konstantin
/// @section LICENSE 
/// This file is part of Onza FDTD.
///
/// Onza FDTD is free software: you can redistribute it and/or modify
/// it under the terms of the GNU General Public License as published by
/// the Free Software Foundation, either version 3 of the License, or
/// (at your option) any later version.
///
/// Onza FDTD is distributed in the hope that it will be useful,
/// but WITHOUT ANY WARRANTY; without even the implied warranty of
/// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
/// GNU General Public License for more details.
///
/// You should have received a copy of the GNU General Public License
/// along with Onza FDTD.  If not, see <http://www.gnu.org/licenses/>.
/// @date   Wed Apr 25 15:01:02 2012
///
/// @brief  Top level calls for FDTD simulation, Doxygen mainpage description.
#include <mpi.h>
#include <cstdio>
#include "mpi-decomposition/halo-exchange-process.h"
/// @brief Init, run and complete simulation.
///
/// @param argc
/// @param argv well known input parameters.
///
/// @return Zero by default.
int main(int argc, char *argv[]) {
  MPI_Init(&argc, &argv);
  onza::HaloExchangeProcess halo_exchange_process;
  int done_status = onza::kDone;
  while (1) {  // use break to report error with done_status.
    done_status = halo_exchange_process.Init(argc, argv);
    if (done_status != onza::kDone) break;
    // Split total simualation domain into subdomains for parallel
    // execution. Each subdomain is simulated on its own MPI processes
    // (physicaly it would be a single core, processor, node etc.).
    done_status = halo_exchange_process.RunDecomposition();
    if (done_status != onza::kDone) break;
    done_status = halo_exchange_process.InitSimulation();
    if (done_status != onza::kDone) break;
    done_status = halo_exchange_process.RunSimulation();
    break;
  }  // end of while breaked with errors
  MPI_Finalize();
  return done_status;
}  // end of main
// ************************************************************************* //
/// \page DevPlan Plan of the Development (in Russian)
/// \tableofcontents
/// Описание необходимых возможностей программы, они же - черновые
/// планы для работы.
/// \section total Предполагаемая схема работы программы в целом.
///
/// Возможно, когда-то к программе будет графическая оболочка. Тогда
/// создание всех управлящих файлов и последовательный вызов
/// подпрограмм будет на её совести. До этого - наборы скриптов и
/// ручная правка конфигов.
/// -# Пишем входной конфиг - в нем задается геометрия модели,
/// положение источников, точки мониторы (в которых снимамются
/// значения), поверхности (через которые считаются, например,
/// потоки), материалы (наример, дисперсионные соотношения для
/// металлов) и т.д.
/// -# Запускаем onza с этим конфигом.
/// -# Получаем на выходе набор файлов, который запрашивали в конфиге.
/// -# Какая-то пост-обработка результатов (рисуем графики, картинки,
/// видео).
/// \section start Ближайшие планы
/// * Объяснить dmitry, что делает программа в целом.
/// * Добавить скрипты, для hg (push, pull)
/// * Просмотреть код, найти места, которые не очень понятны. Либо
///  поменять код, либо добавить комментариев.
/// * \ref selftests
/// * \ref Boundary
/// * \ref geometry
/// * \ref config
/// * Вычистить код, что бы тот проходил cpplint.
/// * Сделать ведение логов и обработчик ошибок.
/// * Вычистить //debug
///
/// \subsection config Конфигуарционный файл.
///
/// Все настройки, которые необходимы для работы программы в едином
/// файле, возможно с ХМL разметкой, а может что-то более
/// простое. Желательно, чтобы можно было решать задачи оптимизации по
/// одному или нескольким параметрам. Или запускать сканирование по
/// меняющемуся параметру (в диапазоне таком-то сделать только-то
/// прогонов, либо в диапазоне таком-то сканировать с таким шагом,
/// либо сканировать по параметру начиная с этого значения с таким-то
/// шагом столько-то раз).
///
/// \subsection geometry Задание геометрии модели.
///
/// Задание геометрии - аналитическое. Т.е. "Хочу круг радиуса
/// столько-то там-то, с таким-то распределением eps от угла поворота
/// и радиуса" + положение и форма источников + мониторы. Нужна некая
/// совместимость с пакетами 3D моделирования, которые понимают
/// скрипты. Например, подпрограммка, которая из конфига для onza
/// делает скрипт на питоне, а этот скрипт можно запустить в FreeCAD,
/// где и будет наглядно показана геометрия.
///
/// \subsection output Вывод результатов.
///
/// Вывод из программы: через файлы формата HDF5. Это касается как и
/// полных snapshot-ов полей, так и спектров, потоков и пр. Для полей
/// - должна быть опция прореженного вывода (например, хочу вывоводить
/// напряжённость электрического поля, 20 точек на микрон).
///
/// \subsection selftests Набор тестовых задач.
///
///  Для автоматической проверки корректности работы программы.
/// Нужны простые задачи, что бы было понятно, не поломалось ли чего
/// из-за внесенных изменений. Другими словами, перед каждым push на
/// общий сервер скрипт, который делает push, запускает
/// ./selftests.sh, если есть какие-то ошибки, то push не выполняется.
///
/// Примеры тестов:
///
/// * 1D модель, фиксированные параметры числа шагов, размера сетки,
/// вывода. Проверяем, что в выводе какое-то число (например поле в
/// третьей точке) совпадает с контрольным.
///
/// * Когда будет считала потоков - можно сделать проверку уравнения
/// Максвелла (поток через через замкнутую поверхность должен быть
/// равен нулю, когда внутри нет источников\поглащения).
///
/// \subsection Boundary Граничные условия.
///
/// Сделать граничные условия. Т.к. оптимизация с умножением на 0 не
/// работает, то нужны отдельные RunAlgorithm для каждого случая +
/// вначале надо выяснять, содержит ли subdomain текущего процесса
/// область, где надо эти уравнения применять. В идеале, хотелось бы
/// возможности наложения заранее заданных граничных условий (вроде
/// PML) на любые конечно-разностные уравнения (не известно заранее
/// какие).
///
/// Возможно в идеале в качестве поглощаеющего условия так будет
/// работать только CPML.
///
/// В целом должна быть функциональность:
/// * Периодические гран. условия.
/// * PML
/// * АBC
/// * UPML
/// * CPML
/// * PEC
/// * PMC
///
/// \subsection Metal Металлы.
///
///  Для начала Друде, Лоренц. Потом более сложные модели, типа PLRC,
///  может можно задавать что-то приближенное к реальности.
///
/// \section Future Когда-нибудь
///
/// Например, \ref FutureTips
// /// This page contains the subsections \ref subsection1 and \ref subsection2.
// /// For more info see page \ref page2.
///
// ************************************************************************* //
/// @page FutureTips Tips for future
/// * To deal with staircase use coordinate transform from:
///  - Title: A General Procedure for Introducing Structured Nonorthogonal
/// Discretization Grids Into High-Order Finite-Difference Time-Domain Methods
///  - Author(s): Armenta, Roberto B.;Sarris, Costas D.
///  - Source: IEEE TRANSACTIONS ON MICROWAVE THEORY AND TECHNIQUES
///  - Volume: 58  Issue: 7  Pages: 1818-1829  DOI: 10.1109/TMTT.2010.2049921
// ************************************************************************* //
/// @page TopLevelAlgorithm Top level algorithm steps
/// - Start exchange process (onza::HaloExchangeProcess::Init()). For each
///   MPI process:
///   * Get MPI runtime parameters
///   * Read simulation global input config
///     (onza::SimulationInputConfig::ReadConfig()).
///   * Cary out domain decomposition
///     onza::HaloExchangeProcess::RunDecomposition().  For early stages
///     of development we will assume star topology, decomposition is
///     optimized do minize volume of exchange data.
/// - Start simulation on subdomain with help of onza::BasicSimulationCore
///   * Note: division could be several times longer than
///     multiplication. This way do not use division while time
///     stepping. All divisions should be done at inititalizing step.
/// @todo2 (For all project) Should be some logging system in Onza.
/// @todo3 (For all project) Perform out of range check for all ceil and float
/// operations.
// ************************************************************************* //
/// @page ChangeLog
/// #ChangeLog
/// ## Version 0.0.2
/// - Implemented 1D, 2D and 3D stepping algorithms.
/// - Some simple input (hardcoded at compilation time).
/// - Adding doxygen related pages.
/// - Speedup at Infiniband cluster 2D checked - 24 times at 16 nodes!
///
/// ## Version 0.0.1
/// - CMake files are configured for ease of use.
// ************************************************************************* //
/// @page "Coding Style"
/// [Google style]: http://google-styleguide.googlecode.com/svn/trunk/cppguide.xml
/// #Coding Style
/// @section Main
/// Use [Google style]!
///
/// @section Additional
/// - Allways end class (methods, functions, namespaces, etc.) declaration
///   with comment
///
///        class GridInputConfig {
///          ...
///        };  // end of class GridInputConfig
///        int  BasicSimulationCore::Init() {
///          ...
///        }  // end of BasicSimulationCore::Init
/// - All classes, containing input parameters for simulation should end
///   with *InputConfig*
// ************************************************************************* //
// /// \page page2 New page
// /// \tableofcontents
// /// Leading text.
// /// \section sec An example section
// /// This page contains the subsections \ref subsection1 and \ref subsection2.
// /// For more info see page \ref page2.
// /// \subsection subsection1 The first subsection
// /// Text.
// /// \subsection subsection2 The second subsection
// /// More text.
// ///
// ************************************************************************* //
/// @namespace onza
/// @brief Top level namespace of Onza FDTD project.
///
/// Also using names from MPI, blitz and std namespaces.
// ************************************************************************* //
// ************************************************************************* //
// ************************************************************************* //
/// @mainpage Onza FDTD Documentation
///
/// [FDTD]: http://en.wikipedia.org/wiki/Finite-difference_time-domain_method
///         "Go to Wikipedia"
/// [Onza source]: http://onzafdtd.org/onza-daily.tar.gz "source"
/// [OpenMPI]: http://www.open-mpi.org
/// [Metamaterials Laboratory]: http://phoi.ifmo.ru/metamaterials
/// [NRU ITMO]: http://en.ifmo.ru "National Research University ITMO"
/// [Blitz++]: http://www.oonumerics.org/blitz
/// [CMake]: http://www.cmake.org "CMake"
/// [Ioffe Institute]: http://www.ioffe.ru/index_en.html
///                    "Ioffe Physical Technical Instute"
///
/// #Parallel generic electromagnetic simulation software
///
/// Onza FDTD is a free (GPLv3+, [Onza source]) high performance
/// electromagnetic simulation software using [finite-difference
/// time-domain (FDTD) method][FDTD]. FDTD method generates very high
/// computational load. To get simulation results in reasonable time
/// your should use a powerful computer or a supercomputer
/// cluster. Onza FDTD was designed to run efficiently with: -
/// multi-core processors - multiprocessor systems - clusters and
/// supercomputers
///
/// Onza FDTD software is developed by [Metamaterials Laboratory] of
/// Photonics and Optical Informatics Department of [NRU ITMO]. It is
/// also supported by [Ioffe Institute].
///
/// Contact Ladutenko Konstantin <kostyfisik at gmail (.)  com>
/// with any questions about Onza FDTD.
///
/// #Prerequisites
/// - **MPI** (to use Onza FDTD, to compile it) <br><br>
///   Onza FDTD needs some MPI realization to be compiled and to be
///   run. It was developed using [OpenMPI]. For Debian/Ubuntu systems
///   it can be installed with:
///
///         # apt-get install openmpi-bin openmpi-doc libopenmpi-dev
///
///   <br>
/// - **The [Blitz++] library** (to compile Onza FDTD) <br><br>
///   For Debian/Ubuntu systems it can be installed with
///
///         # apt-get install libblitz0-dev
///
///   <br>
/// - **And [CMake] build system** (to compile Onza FDTD) <br><br>
///   Use it on Linux/Unix and other operation systems.
