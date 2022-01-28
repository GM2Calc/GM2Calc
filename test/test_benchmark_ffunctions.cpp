#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN 1

#include "doctest.h"
#include "gm2_ffunctions.hpp"
#include "stopwatch.hpp"

#include <random>

namespace {

/// random number generator
double random(double start, double stop)
{
   static std::minstd_rand gen;
   std::uniform_real_distribution<double> dist(start, stop);
   return dist(gen);
}

/// measure time for calling the function f()
template <class F>
double time_in_milliseconds(F&& f)
{
   gm2calc::Stopwatch sw;
   sw.start();
   f();
   sw.stop();
   return sw.get_time_in_milliseconds();
}

template <class F>
double bench_1_in_ms(double start, double stop, unsigned N, F f)
{
   std::vector<double> x(N);

   const auto ran = [&] { return random(start, stop); };

   std::generate(std::begin(x), std::end(x), ran);

   const auto time_in_ms = time_in_milliseconds(
      [&] {
         for (unsigned i = 0; i < N; ++i) {
            (void) f(x[i]);
         }
      }
   );

   return time_in_ms/N;
}

} // anonymous namespace

TEST_CASE("benchmark f_PS")
{
   const auto time_in_ms = bench_1_in_ms(
      0.0, 0.5, 1000000,
      [] (double x) { return gm2calc::f_PS(x); });

   std::cout << "f_PS(x): average time per point: " << time_in_ms*1000 << " ns\n";
}

TEST_CASE("benchmark F1")
{
   const auto time_in_ms = bench_1_in_ms(
      0.1, 0.5, 1000000,
      [] (double x) { return gm2calc::F1(x); });

   std::cout << "F1(x): average time per point: " << time_in_ms*1000 << " ns\n";
}

TEST_CASE("benchmark F1t")
{
   const auto time_in_ms = bench_1_in_ms(
      0.1, 0.5, 1000000,
      [] (double x) { return gm2calc::F1t(x); });

   std::cout << "F1t(x): average time per point: " << time_in_ms*1000 << " ns\n";
}

TEST_CASE("benchmark F2")
{
   const auto time_in_ms = bench_1_in_ms(
      0.1, 0.5, 1000000,
      [] (double x) { return gm2calc::F2(x); });

   std::cout << "F2(x): average time per point: " << time_in_ms*1000 << " ns\n";
}

TEST_CASE("benchmark F3")
{
   const auto time_in_ms = bench_1_in_ms(
      0.1, 0.5, 1000000,
      [] (double x) { return gm2calc::F3(x); });

   std::cout << "F3(x): average time per point: " << time_in_ms*1000 << " ns\n";
}

TEST_CASE("benchmark Iabc")
{
   const unsigned N = 1000000;
   std::vector<double> x(N), y(N) ,z(N);

   const auto ran = [] { return random(0.1, 1000); };

   std::generate(std::begin(x), std::end(x), ran);
   std::generate(std::begin(y), std::end(y), ran);
   std::generate(std::begin(z), std::end(z), ran);

   const auto time_in_ms = time_in_milliseconds(
      [&] {
         for (unsigned i = 0; i < N; ++i) {
            (void) gm2calc::Iabc(x[i], y[i], z[i]);
         }
      }
   );

   std::cout << "Iabc(x,y,z): average time per point: " << time_in_ms*1000/N << " ns\n";
}

TEST_CASE("benchmark Phi")
{
   const unsigned N = 1000000;
   std::vector<double> x(N), y(N) ,z(N);

   const auto ran = [] { return random(0.1, 1000); };

   std::generate(std::begin(x), std::end(x), ran);
   std::generate(std::begin(y), std::end(y), ran);
   std::generate(std::begin(z), std::end(z), ran);

   const auto time_in_ms = time_in_milliseconds(
      [&] {
         for (unsigned i = 0; i < N; ++i) {
            (void) gm2calc::Phi(x[i], y[i], z[i]);
         }
      }
   );

   std::cout << "Phi(x,y,z): average time per point: " << time_in_ms*1000/N << " ns\n";
}
