// [[Rcpp::plugins("cpp11")]]
// [[Rcpp::depends(RcppArmadillo, RcppProgress)]]
#include "utils.h"
#include "vips/vips8"
#include <RcppArmadillo.h>
#include <cassert>
#include <chrono>
#include <cmath>
#include <csignal>
#include <filesystem>
#include <indicators/cursor_control.hpp>
#include <indicators/progress_bar.hpp>
#include <string>
#include <tuple>

#ifdef _OPENMP
#include <omp.h>
#endif

bool
is_in_circle(
    const double point_x, const double point_y, const double center_x,
    const double center_y, const double radius
) {
  const double dist = std::sqrt(
      std::pow(point_x - center_x, 2) + std::pow(point_y - center_y, 2)
  );
  return dist <= radius;
}

// The origin point of the coordinate system is at the center of the circle.
bool
is_in_fan(
    const double point_x, const double point_y, const double center_x,
    const double center_y, const double fan_point_1_x,
    const double fan_point_1_y, const double fan_point_2_x,
    const double fan_point_2_y, const double radius
) {
  // The first point of the subspot is left.
  const double diff_1_x =
      (point_y - fan_point_1_y) * (fan_point_1_x - fan_point_2_x) +
      (fan_point_1_x - point_x) * (fan_point_1_y - fan_point_2_y);
  const double diff_1_y =
      (center_y - fan_point_1_y) * (fan_point_1_x - fan_point_2_x) +
      (fan_point_1_x - center_x) * (fan_point_1_y - fan_point_2_y);

  if (diff_1_x * diff_1_y < 0 &&
      !is_in_circle(point_x, point_y, center_x, center_y, radius))
    return false;

  // The second point of the subspot is left.
  const double diff_2_x = (point_y - center_y) * (fan_point_2_x - center_x) +
                          (center_x - point_x) * (fan_point_2_y - center_y);
  const double diff_2_y =
      (fan_point_1_y - center_y) * (fan_point_2_x - center_x) +
      (center_x - fan_point_1_x) * (fan_point_2_y - center_y);

  if (diff_2_x * diff_2_y < 0)
    return false;

  // The third point of the subspot is left.
  const double diff_3_x = (point_y - center_y) * (fan_point_1_x - center_x) +
                          (center_x - point_x) * (fan_point_1_y - center_y);
  const double diff_3_y =
      (fan_point_2_y - center_y) * (fan_point_1_x - center_x) +
      (center_x - fan_point_2_x) * (fan_point_1_y - center_y);

  return diff_3_x * diff_3_y >= 0;
}

bool
get_tile_paths(
    const std::string &barcode, const std::filesystem::path &spot_output_path,
    const std::filesystem::path &subspot_output_path,
    std::vector<std::tuple<std::string, std::filesystem::path>> &tile_paths
) {
  tile_paths[0] = std::make_pair(
      barcode, spot_output_path / std::filesystem::path(barcode + ".png")
  );

  for (int i = 1; i < 7; i++)
    tile_paths[i] = std::make_pair(
        barcode + ":" + std::to_string(i),
        subspot_output_path /
            std::filesystem::path(barcode + "_" + std::to_string(i) + ".png")
    );

  for (auto &i : tile_paths)
    if (!is_regular_file(std::get<1>(i)))
      return false;

  return true;
}

// 1. Rotate other fans to the position of fan_1.
// 2. Rotation by 60 and 120 degrees will scale the output image up by
// (1+sqrt(3))/2 to keep the size of original image unchanged.
// 3. Return: mask image object, clockwise rotation angle, top-left corner
std::vector<std::tuple<vips::VImage, double, int>>
prepare_masks(const int width, const int num_bands = 3) {
  vips::VImage ori = vips::VImage::black(
      width, width, vips::VImage::option()->set("bands", num_bands)
  );
  const int rotated_corner = static_cast<int>((std::sqrt(3) + 1) * width / 4);

  const std::vector<double> white_pxl(num_bands, 255);

  vips::VImage circle, fan_1, fan_2, fan_3, fan_4, fan_5, fan_6;
  circle = ori.copy_memory();
  fan_1  = ori.copy_memory();
  fan_2  = ori.copy_memory();
  fan_3  = ori.copy_memory();
  fan_4  = ori.copy_memory();
  fan_5  = ori.copy_memory();
  fan_6  = ori.copy_memory();

  for (int col = 0; col < width; col++) {
    for (int row = 0; row < width; row++) {
      if (is_in_circle(col, row, width / 2, width / 2, width / 2))
        circle.draw_rect(white_pxl, col, row, 1, 1);

      if (is_in_fan(
              col, row, width / 2, width / 2, width / 2, width,
              (2 + std::sqrt(3)) * width / 4, 3 * width / 4, width / 2
          ))
        fan_1.draw_rect(white_pxl, col, row, 1, 1);

      if (is_in_fan(
              col, row, width / 2, width / 2, width / 2, width,
              (2 - std::sqrt(3)) * width / 4, 3 * width / 4, width / 2
          ))
        fan_2.draw_rect(white_pxl, col, row, 1, 1);

      if (is_in_fan(
              col, row, width / 2, width / 2, width / 2, 0,
              (2 + std::sqrt(3)) * width / 4, width / 4, width / 2
          ))
        fan_3.draw_rect(white_pxl, col, row, 1, 1);

      if (is_in_fan(
              col, row, width / 2, width / 2, width / 2, 0,
              (2 - std::sqrt(3)) * width / 4, width / 4, width / 2
          ))
        fan_4.draw_rect(white_pxl, col, row, 1, 1);

      if (is_in_fan(
              col, row, width / 2, width / 2, (2 + std::sqrt(3)) * width / 4,
              width / 4, (2 + std::sqrt(3)) * width / 4, 3 * width / 4,
              width / 2
          ))
        fan_5.draw_rect(white_pxl, col, row, 1, 1);

      if (is_in_fan(
              col, row, width / 2, width / 2, (2 - std::sqrt(3)) * width / 4,
              width / 4, (2 - std::sqrt(3)) * width / 4, 3 * width / 4,
              width / 2
          ))
        fan_6.draw_rect(white_pxl, col, row, 1, 1);
    }
  }

  return std::vector<std::tuple<vips::VImage, double, int>>{
      std::tuple<vips::VImage, double, int>{circle, 0, 0},
      std::tuple<vips::VImage, double, int>{
          fan_1, 0, static_cast<int>(width / 2)
      },
      std::tuple<vips::VImage, double, int>{fan_2, -60, rotated_corner},
      std::tuple<vips::VImage, double, int>{fan_3, 120, rotated_corner},
      std::tuple<vips::VImage, double, int>{
          fan_4, 180, static_cast<int>(width / 2)
      },
      std::tuple<vips::VImage, double, int>{fan_5, 60, rotated_corner},
      std::tuple<vips::VImage, double, int>{fan_6, -120, rotated_corner}
  };
}

std::vector<std::tuple<int, int>>
slice_data(const int num_to_slice, const int num_chunks) {
  if (num_to_slice <= num_chunks)
    return std::vector<std::tuple<int, int>>{std::make_pair(0, num_chunks - 1)};

  std::vector<std::tuple<int, int>> ret;

  const int len = num_to_slice / num_chunks;
  int mod       = num_to_slice % num_chunks;
  int start = 0, end = len - 1;

  for (int i = 0; i < num_chunks; i++) {
    if (mod > 0) {
      mod--;
      end++;
    }

    ret.emplace_back(std::make_pair(start, end));

    start = end + 1;
    end   = start + len - 1;
  }

  return ret;
}

template <typename T>
void
flatten_image(
    const vips::VImage &img, arma::Mat<T> &flattened_img, const size_t col_idx,
    std::mutex *__mutex = nullptr
) {
  size_t buf_len = 0;

  // Must always have an unsigned char pointer returned.
  // The image is flattened in such an order: top to bottom, left to right,
  // bands.
  const unsigned char *flattened =
      (unsigned char *) img.write_to_memory(&buf_len);

  std::vector<uint32_t> reordered(buf_len);

  // The flattened image should be reordered in: bands, left to right, top
  // to bottom
  size_t i = 0;
  for (int __band = 0; __band < img.bands(); __band++)
    for (int __row = 0; __row < img.height(); __row++)
      for (int __col = 0; __col < img.width(); __col++)
        reordered[i++] =
            static_cast<uint32_t>(flattened
                                      [__row * img.width() * img.bands() +
                                       __col * img.bands() + __band]);

  {
    if (__mutex != nullptr)
      std::lock_guard<std::mutex> lock(*__mutex);

    flattened_img.col(col_idx) = arma::Col<T>(reordered);
  }

  g_free((void *) flattened);
  flattened = nullptr;
}

void
__get_spot_subspot_tiles_from_image(
    const Rcpp::CharacterVector &barcodes,
    const arma::mat &spot_center_coordinates, const double spot_radius_pxl,
    const std::string &fullres_image_file, const std::string &tile_image_dir,
    arma::umat &spot_flatten_images, arma::umat &subspot_flatten_images,
    std::vector<std::string> &subspot_barcodes, const int thread_num,
    const bool verbose = false
) {
  std::vector<int> thread_hits;

#ifdef _OPENMP
  omp_set_max_active_levels(2);
  omp_set_num_threads(thread_num);

  for (int i = 0; i < thread_num; i++)
    thread_hits.emplace_back(0);

  if (verbose) {
    std::cout << "[DEBUG] The number of threads is " << thread_num << std::endl;
  }
#endif

  const std::filesystem::path __spot_output_path(
      tile_image_dir / std::filesystem::path("spot")
  );
  std::filesystem::create_directories(__spot_output_path);

  const std::filesystem::path __subspot_output_path(
      tile_image_dir / std::filesystem::path("subspot")
  );
  std::filesystem::create_directories(__subspot_output_path);

  const int __spot_radius_pxl = std::floor(spot_radius_pxl);

  // Load image.
  const vips::VImage img =
      vips::VImage::new_from_file(fullres_image_file.c_str());

  // Prepare for the outputs.
  spot_flatten_images = arma::umat(
      img.bands() * std::pow(2 * __spot_radius_pxl, 2), barcodes.length()
  );
  subspot_flatten_images = arma::umat(
      img.bands() * std::pow(__spot_radius_pxl, 2), 6 * barcodes.length()
  );
  subspot_barcodes = std::vector<std::string>(6 * barcodes.length());

  // Path to spot and subspot tiles.
  std::vector<std::tuple<std::string, std::filesystem::path>> tile_paths(7);

  // Masks for spot and subspot images.
  const std::vector<std::tuple<vips::VImage, double, int>> masks =
      prepare_masks(2 * __spot_radius_pxl, img.bands());

  // Progree bar.
  indicators::show_console_cursor(false);
  indicators::ProgressBar pb{
      indicators::option::MaxProgress{barcodes.length() - 1},
      indicators::option::BarWidth{50},
      indicators::option::Start{" ["},
      indicators::option::Fill{"█"},
      indicators::option::Lead{"█"},
      indicators::option::Remainder{"-"},
      indicators::option::End{"]"},
      indicators::option::PrefixText{"Slicing"},
      indicators::option::ForegroundColor{indicators::Color::blue},
      indicators::option::ShowElapsedTime{true},
      indicators::option::ShowRemainingTime{true},
      indicators::option::FontStyles{
          std::vector<indicators::FontStyle>{indicators::FontStyle::bold}
      }
  };

#pragma omp parallel for
  for (int i = 0; i < barcodes.length(); i++) {
    pb.tick();

#ifdef _OPENMP
#pragma omp atomic update
    thread_hits[omp_get_thread_num()]++;
#endif

    if (get_tile_paths(
            static_cast<std::string>(barcodes[i]), __spot_output_path,
            __subspot_output_path, tile_paths
        )) {
      for (int j = 0; j < 7; j++) {
        if (j > 0)
          subspot_barcodes[(j - 1) * barcodes.length() + i] =
              std::get<0>(tile_paths[j]);

        flatten_image(
            vips::VImage::new_from_file(std::get<1>(tile_paths[j]).c_str()),
            j == 0 ? spot_flatten_images : subspot_flatten_images,
            j == 0 ? i : (j - 1) * barcodes.length() + i
        );
      }
    } else {
      const int center_col = spot_center_coordinates(i, 0);
      const int center_row = spot_center_coordinates(i, 1);

      const int left_col = center_col - __spot_radius_pxl;
      const int top_row  = center_row - __spot_radius_pxl;

      const vips::VImage __rect = img.extract_area(
          left_col, top_row, 2 * __spot_radius_pxl, 2 * __spot_radius_pxl
      );

      for (int j = 0; j < 7; j++) {
        vips::VImage __masked =
            (__rect & std::get<0>(masks[j]))
                .similarity(
                    vips::VImage::option()->set("angle", std::get<1>(masks[j]))
                );

        if (j > 0) {
          __masked = __masked.extract_area(
              std::get<2>(masks[j]), std::get<2>(masks[j]), __spot_radius_pxl,
              __spot_radius_pxl
          );

          subspot_barcodes[(j - 1) * barcodes.length() + i] =
              std::get<0>(tile_paths[j]);
        }

        flatten_image(
            __masked, j == 0 ? spot_flatten_images : subspot_flatten_images,
            j == 0 ? i : (j - 1) * barcodes.length() + i
        );

        __masked.write_to_file(std::get<1>(tile_paths[j]).c_str());
      }
    }
  }

  indicators::show_console_cursor(true);

#ifdef _OPENMP
  if (verbose) {
    print_thread_hits(thread_hits);
  }
#endif
}

// [[Rcpp::export]]
Rcpp::List
get_spot_subspot_tiles_from_image(
    const Rcpp::CharacterVector &barcodes,
    const arma::mat &spot_center_coordinates, const double spot_radius_pxl,
    const std::string &fullres_image_file, const std::string &tile_image_dir,
    const bool init_vips = true, const bool shutdown_vips = true,
    const int thread_num = 1, const bool verbose = false
) {
  assertm(
      barcodes.length() == spot_center_coordinates.n_rows,
      "The number of spots in `barcodes` and `spot_center_coordinates` do "
      "not match."
  );

  arma::umat spot_flatten_images;
  arma::umat subspot_flatten_images;
  std::vector<std::string> subspot_barcodes;

  if (init_vips)
    vips_init("");

  __get_spot_subspot_tiles_from_image(
      barcodes, spot_center_coordinates, spot_radius_pxl, fullres_image_file,
      tile_image_dir, spot_flatten_images, subspot_flatten_images,
      subspot_barcodes, thread_num, verbose
  );

  if (shutdown_vips)
    vips_shutdown();

  return Rcpp::List::create(
      Rcpp::_["spot"]             = spot_flatten_images,
      Rcpp::_["subspot"]          = subspot_flatten_images,
      Rcpp::_["subspot_barcodes"] = subspot_barcodes
  );
}
