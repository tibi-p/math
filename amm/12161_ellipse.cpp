#include <iostream>
#include <vector>

void compute_stats(int limit)
{
  std::vector<std::vector<std::pair<int, int>>> sums(limit + 1);
  for (int i = 0; ; i++) {
    int partial = i * i;
    if (partial > limit) {
      break;
    }
    for (int j = 0; ; j++) {
      int result = partial + 3 * j * j;
      if (result > limit) {
        break;
      }
      sums[result].push_back(std::make_pair(i, j));
    }
  }
  for (int i = 0; i <= limit; i += 4) {
    if (sums[i].empty()) {
      continue;
    }
    std::cout << i << " (" << sums[i].size() << "): ";
    for (auto elem: sums[i]) {
      std::cout << " (" << elem.first << ',' << elem.second << ')';
    }
    std::cout << std::endl;
  }
}

int main()
{
  compute_stats(10000);
  return 0;
}
