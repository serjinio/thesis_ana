#include <iostream>
#include <vector>
#include <algorithm>


int bad_neighbors_subsqeq(std::vector<int> input) {
  std::vector<int> partial_sums;
  for (size_t i = 0; i < input.size(); ++i) {
    partial_sums.push_back(input[i]);
  }

  for (size_t i = 1; i < input.size(); ++i) {
    int max_donation = 0;
    for (size_t k = 0; k <= i - 1; ++k) {
      if (i - k == 1) {
        continue;
      }

      int donation = partial_sums[k] + input[i];
      if (donation > max_donation) {
        max_donation = donation;
        partial_sums[i] = max_donation;
      }
    }
  }

  std::cout << "Input        : ";
  for (size_t i = 0; i < input.size(); ++i) {
    std::cout << input[i] << "; ";
  }
  std::cout << std::endl;
  std::cout << "Donation sums: ";
  for (size_t i = 0; i < input.size(); ++i) {
    std::cout << partial_sums[i] << "; ";
  }
  std::cout << std::endl;

  return *std::max_element(std::begin(partial_sums), std::end(partial_sums));
}

int bad_neighbors(std::vector<int> input) {
  if (input.size() <= 3) {
    return *std::max_element(input.begin(), input.end());
  }

  std::vector<int> subseq1(input.size());
  std::vector<int> subseq2(input.size());
  std::copy(input.begin(), input.end() - 1, subseq1.begin());
  std::copy(input.begin() + 1, input.end(), subseq2.begin());
  return std::max(bad_neighbors_subsqeq(subseq1), bad_neighbors_subsqeq(subseq2));
}

std::vector<int> flower_garden(std::vector<int> height, std::vector<int> bloom,
                               std::vector<int> wilt) {
  std::vector<int> garden_idxes;
  std::vector<int> sorted_indices(height.size());
  for (size_t i = 0; i < height.size(); i++) {
    sorted_indices[i] = i;
  }
  std::sort(sorted_indices.begin(), sorted_indices.end(),
            [&height](int a, int b) -> bool { return height[a] > height[b]; });
  std::cout << "sorted indices: ";
  for (size_t i = 0; i < height.size(); i++) {
    std::cout << sorted_indices[i] << "; ";
  }
  std::cout << std::endl;

  bool blocks[height.size()][height.size()];
  for (size_t i = 0; i < height.size() - 1; ++i) {
    int bi = bloom[sorted_indices[i]];
    int wi = wilt[sorted_indices[i]];
    for (size_t k = 0; k < height.size(); ++k) {
      if (i == k) {
        blocks[i][k] = false;
        continue;
      }
      int bk = bloom[sorted_indices[k]];
      int wk = wilt[sorted_indices[k]];
      if ((bk <= wi || wk >= bi) && i < k) {
        blocks[i][k] = true;
      } else {
        blocks[i][k] = false;
      }
    }
  }
  for (size_t i = 0; i < height.size(); ++i) blocks[height.size() - 1][i] = false;
  std::cout << "blocks array: " << std::endl;
  for (size_t i = 0; i < height.size(); i++) {
    for (size_t k = 0; k < height.size(); k++) {
      std::cout << (blocks[i][k] == true ? 1 : 0) << "; ";
    }
    std::cout << std::endl;
  }

  garden_idxes.push_back(0);
  for (size_t i = 1; i < height.size(); ++i) {
    int new_flower = i;
    int insert_idx = garden_idxes.size();
    for (int k = garden_idxes.size() - 1; k >= 0; --k) {
      std::cout << "block[" << garden_idxes[k] << "][" << i << "]"
                << "==" << blocks[garden_idxes[k]][i] << std::endl;
      if (blocks[garden_idxes[k]][i]) {
        insert_idx = k;
      }
    }
    std::cout << "inserting " << height[sorted_indices[i]]
              << " at index: " << insert_idx << std::endl;
    garden_idxes.insert(garden_idxes.begin() + insert_idx, new_flower);
  }

  std::vector<int> garden_heights(height.size());
  for (size_t i = 0; i < garden_heights.size(); i++) {
    garden_heights[i] = height[garden_idxes[i]];
  }

  return garden_heights;
}

int main()
{
  //std::vector<int> input{ 10, 3, 2, 5, 7, 8 };
  std::vector<int> input{ 94, 40, 49, 65, 21, 21, 106, 80, 92, 81, 679, 4, 61,
      6, 237, 12, 72, 74, 29, 95, 265, 35, 47, 1, 61, 397,
      52, 72, 37, 51, 1, 81, 45, 435, 7, 36, 57, 86, 81, 72 };
  // std::vector<int> input{ 1, 2, 3, 4, 5, 1, 2, 3, 4, 5 };

  std::vector<int> height = {5,4,3,2,1};
  std::vector<int> bloom = {1,1,1,1,1};
  std::vector<int> wilt = {365,365,365,365,365};

  auto res = flower_garden(height, bloom, wilt);
  std::cout << "Flower garden_idxes: ";
  for (size_t i = 0; i < res.size(); i++) {
    std::cout << res[i] << "; ";
  }
  std::cout << std::endl;
}
