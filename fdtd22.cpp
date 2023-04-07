#include <array>
#include <vector>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <cassert>
#include <filesystem>

struct sound_source {
  sound_source(float x, float y, int t) {
    xpos = x;
    ypos = y;
    tstep = t;
    data.resize(t);
  }
  sound_source(float x, float y, std::vector<float>&& raw_data) {
    xpos = x;
    ypos = y;
    data = std::move(raw_data);
    tstep = data.size();
  }
  
  float xpos, ypos, duration;
  int tstep;
  std::vector<float> data;
}; 

template<int x_size, int y_size>
class fdtd_solver {
  public:
    fdtd_solver(float fmax, float duration, float sound_speed = 314, float raw = 0.001, float kappa = 0.1) {
      dt_ = 1 / (2 * fmax);
      h_ = std::sqrt(2) * sound_speed * dt_;
      duration_ = duration;
      raw_ = raw;
      kappa_ = kappa;
      t_step_ = int(duration / dt_);
      vx_ = std::vector<std::array<float, (x_size + 1) * y_size>>(t_step_);
      vy_ = std::vector<std::array<float, x_size * (y_size + 1)>>(t_step_);
      p_ = std::vector<std::array<float, x_size * y_size>>(t_step_);
    }

    void solve(const sound_source& external) {
      assert(external.tstep == t_step_&&"t step is different.");
      int x_ex = static_cast<int>(external.xpos);
      int y_ex = static_cast<int>(external.ypos);

      p_[0][p_id(x_ex, y_ex)] = external.data[0];

      float v_fac = dt_ / (raw_ * h_);
      float p_fac = dt_ * kappa_ / h_;
      // update 
      for (int t = 1; t < t_step_; t++) {
        // update particle velocity
        for (int i = 0; i < x_size; i++) {
          for (int j = 0; j < y_size; j++) {
            // vx
            if (i < x_size - 1) {
              vx_[t][vx_id(i, j)] = 
                vx_[t-1][vx_id(i, j)] - v_fac * (p_[t-1][p_id(i+1, j)] - p_[t-1][p_id(i, j)]);
            }
            // vy
            if (j < y_size - 1) {
              vy_[t][vy_id(i, j)] = 
                vy_[t-1][vy_id(i, j)] - v_fac * (p_[t-1][p_id(i, j+1)] - p_[t-1][p_id(i, j)]);
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
              vx_[t][vx_id(i, j)] - vx_[t][vx_id(i-1, j)] + 
              vy_[t][vy_id(i, j)] - vy_[t][vy_id(i, j-1)]
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
          ofs << std::fixed << std::setprecision(4) << grid[i] << ',';
        }
        ofs << grid[x_size * y_size - 1] << std::endl;
      }

      ofs.close();
    }

    inline float get_dt() const { return dt_; }
    inline int   get_t_step() const { return t_step_; }

  private:
    inline int p_id(int x, int y) { return x + x_size * y; }
    inline int vx_id(int x, int y) { return x + 1 + (x_size + 1) * y; }
    inline int vy_id(int x, int y) { return x + x_size * (y + 1); }

    std::vector<std::array<float, (x_size + 1) * y_size>> vx_;
    std::vector<std::array<float, x_size * (y_size + 1)>> vy_;
    std::vector<std::array<float, x_size * y_size>>       p_;

    float h_; // grid_size
    float dt_;
    float duration_;
    float raw_;
    float kappa_;
    int   t_step_;
};

std::vector<float> create_sin_signal(float amp, float freq, float duration, float rate) {
  std::vector<float> ret;
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
  float fmax = 100;
  float rate = fmax * 2;
  float duration = 1;
  float sound_speed = 10;
  fdtd_solver<xsize, ysize> solver(fmax, duration, sound_speed);
  auto sin = create_sin_signal(128, 50, duration, rate);
  auto pulse = std::vector<float>(duration * rate, 0);
  pulse[0] = 128;
  sound_source source(xsize / 2, ysize / 2, std::move(sin));
  solver.solve(source);
  solver.save_pressure("fdtd22sin128");
}