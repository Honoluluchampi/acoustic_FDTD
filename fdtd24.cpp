#include <array>
#include <vector>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <cassert>
#include <filesystem>

struct sound_source {
  sound_source(double x, double y, int t) {
    xpos = x;
    ypos = y;
    tstep = t;
    data.resize(t);
  }
  sound_source(double x, double y, std::vector<double>&& raw_data) {
    xpos = x;
    ypos = y;
    data = std::move(raw_data);
    tstep = data.size();
  }
  
  double xpos, ypos, duration;
  int tstep;
  std::vector<double> data;
}; 

template<int x_size, int y_size>
class fdtd_solver {
  public:
    fdtd_solver(double fmax, double duration, double sound_speed = 314, double raw = 0.001, double kappa = 0.1) {
      dt_ = 1 / (2 * fmax);
      h_ = std::sqrt(2) * sound_speed * dt_;
      duration_ = duration;
      raw_ = raw;
      kappa_ = kappa;
      t_step_ = int(duration / dt_);
      vx_ = std::vector<std::array<double, (x_size + 3) * y_size>>(t_step_);
      vy_ = std::vector<std::array<double, x_size * (y_size + 3)>>(t_step_);
      p_ = std::vector<std::array<double, x_size * y_size>>(t_step_);
    }

    void solve(const sound_source& external) {
      assert(external.tstep == t_step_&&"t step is different.");
      int x_ex = static_cast<int>(external.xpos);
      int y_ex = static_cast<int>(external.ypos);

      // initial state
      p_[0][p_id(x_ex, y_ex)] = external.data[0];

      double v_fac = dt_ / (raw_ * h_);
      double p_fac = dt_ * kappa_ / h_;

      // taylor factor
      double c0 = 9.f / 8.f;
      double c1 = 1.f / 24.f;

      // update 
      for (int t = 1; t < t_step_; t++) {
        // update particle velocity
        for (int i = 1; i < x_size; i++) {
          for (int j = 1; j < y_size; j++) {
            // vx
            if (i < x_size - 2) {
              vx_[t][vx_id(i, j)] = 
                vx_[t-1][vx_id(i, j)] - v_fac * (
                  c0 * (p_[t-1][p_id(i+1, j)] - p_[t-1][p_id(i, j)]) - 
                  c1 * (p_[t-1][p_id(i+2, j)] - p_[t-1][p_id(i-1, j)]) 
                );
            }
            // vy
            if (j < y_size - 2) {
              vy_[t][vy_id(i, j)] = 
                vy_[t-1][vy_id(i, j)] - v_fac * (
                  c0 * (p_[t-1][p_id(i, j+1)] - p_[t-1][p_id(i, j)]) - 
                  c1 * (p_[t-1][p_id(i, j+2)] - p_[t-1][p_id(i, j-1)]) 
                );
            }
          }
        }
        // update pressure
        for (int i = 0; i < x_size; i++) {
          for (int j = 0; j < y_size; j++) {
            // add external force
            if (i == x_ex && j == y_ex) {
              p_[t][p_id(i, j)] = external.data[t];
              continue;
            }

            p_[t][p_id(i, j)] = p_[t-1][p_id(i, j)] - p_fac * (
              c0 * (vx_[t][vx_id(i, j)] - vx_[t][vx_id(i-1, j)] + vy_[t][vy_id(i, j)] - vy_[t][vy_id(i, j-1)]) -
              c1 * (vx_[t][vx_id(i+1, j)] - vx_[t][vx_id(i-2, j)] + vy_[t][vy_id(i, j+1)] - vy_[t][vy_id(i, j-2)])
            );
          }
        }
      }
    }

    void save_pressure(std::string name = "output") {
      if (std::filesystem::exists(name))
        std::filesystem::create_directory(name);
      std::ofstream ofs("output.txt");
      // constants
      ofs << name << std::endl;
      ofs << x_size << std::endl;
      ofs << y_size << std::endl;

      for (const auto& grid : p_) {
        for (int i = 0; i < x_size * y_size - 1; i++) {
          ofs << std::fixed << std::setprecision(3) << grid[i] << ',';
        }
        ofs << grid[x_size * y_size - 1] << std::endl;
      }

      ofs.close();
    }

    inline double get_dt() const { return dt_; }
    inline int   get_t_step() const { return t_step_; }

  private:
    inline int p_id(int x, int y) { return x + x_size * y; }
    inline int vx_id(int x, int y) { return x + 2 + (x_size + 3) * y; }
    inline int vy_id(int x, int y) { return x + x_size * (y + 2); }

    std::vector<std::array<double, (x_size + 3) * y_size>> vx_;
    std::vector<std::array<double, x_size * (y_size + 3)>> vy_;
    std::vector<std::array<double, x_size * y_size>>       p_;

    double h_; // grid_size
    double dt_;
    double duration_;
    double raw_;
    double kappa_;
    int   t_step_;
};

std::vector<double> create_sin_signal(double amp, double freq, double duration, double rate) {
  std::vector<double> ret;
  int tstep = rate * duration;
  ret.resize(tstep);
  for (int t = 0; t < tstep; t++) {
    ret[t] = amp * std::sin(2 * M_PI * freq * t / rate);
  }
  return ret;
}

int main() {
  const int xsize = 128;
  const int ysize = 128;
  double fmax = 100;
  double rate = fmax * 2;
  double duration = 1;
  double sound_speed = 10;
  fdtd_solver<xsize, ysize> solver(fmax, duration, sound_speed);
  auto sin = create_sin_signal(128, 50, duration, rate);
  auto pulse = std::vector<double>(duration * rate, 0);
  pulse[0] = 128;
  sound_source source(xsize / 2, ysize / 2, std::move(pulse));
  solver.solve(source);
  solver.save_pressure("fdtd24pulse128");
}