#ifndef SRC_ADVENTURE_H_
#define SRC_ADVENTURE_H_

#include <algorithm>
#include <vector>

#include "../third_party/threadpool/threadpool.h"

#include "./types.h"
#include "./utils.h"

class Adventure {
 public:
  virtual ~Adventure() = default;

  virtual uint64_t packEggs(std::vector<Egg> eggs, BottomlessBag &bag) = 0;

  virtual void arrangeSand(std::vector<GrainOfSand> &grains) = 0;

  virtual Crystal selectBestCrystal(std::vector<Crystal> &crystals) = 0;
};

static std::vector<GrainOfSand>::iterator part(
    std::vector<GrainOfSand>::iterator start,
    std::vector<GrainOfSand>::iterator end) {
  auto res = std::prev(end, 1);

  auto i = start;

  for (auto j = start; j != res; ++j) {
    if (*j < *res) {
      std::swap(*i++, *j);
    }
  }

  std::swap(*i, *res);

  return i;
}

/* Sequential implementation of 3 algorithms. */
class LonesomeAdventure : public Adventure {
 public:
  LonesomeAdventure() = default;

  virtual uint64_t packEggs(std::vector<Egg> eggs, BottomlessBag &bag) {
    /* return packHelper with pos on the last elem in eggs vector. */
    std::vector<std::vector<uint64_t>> mem;

    /* Fill matrix 'mem' with 0 values. */
    for (uint64_t i = 0; i <= eggs.size(); i++) {
      std::vector<uint64_t> v;

      for (uint64_t j = 0; j <= bag.getCapacity(); j++) {
        v.push_back(0);
      }

      mem.push_back(v);
    }

    return packHelper(eggs, bag.getCapacity(), eggs.size(), mem);
  }

  virtual void arrangeSand(std::vector<GrainOfSand> &grains) {
    quickSort(grains.begin(), grains.end());
  }

  virtual Crystal selectBestCrystal(std::vector<Crystal> &crystals) {
    Crystal maxCrystal(0);

    for (const auto &crystal : crystals) {
      if (maxCrystal < crystal) maxCrystal = crystal;
    }

    return maxCrystal;
  }

 private:
  /* Helper function for egg packing with dynamic programming. */
  static uint64_t packHelper(std::vector<Egg> &eggs, uint64_t capacity,
                             uint64_t pos,
                             std::vector<std::vector<uint64_t>> &mem) {
    if (pos == 0) {
      mem[pos][capacity] = 0;

      return mem[pos][capacity];
    }

    if (mem[pos][capacity] != 0) return mem[pos][capacity];

    if (eggs[pos - 1].getSize() > capacity) {
      mem[pos][capacity] = packHelper(eggs, capacity, pos - 1, mem);

      return mem[pos][capacity];
    }

    mem[pos][capacity] = std::max(
        eggs[pos - 1].getWeight() +
            packHelper(eggs, capacity - eggs[pos - 1].getSize(), pos - 1, mem),
        packHelper(eggs, capacity, pos - 1, mem));

    return mem[pos][capacity];
  }

  static void quickSort(std::vector<GrainOfSand>::iterator start,
                        std::vector<GrainOfSand>::iterator end) {
    if (std::distance(start, end) > 1) {
      auto q = part(start, end);
      quickSort(start, q);
      quickSort(q + 1, end);
    }
  }
};

#define THRESHOLD 200

class TeamAdventure : public Adventure {
 public:
  explicit TeamAdventure(uint64_t numberOfShamansArg)
      : numberOfShamans(numberOfShamansArg),
        councilOfShamans(numberOfShamansArg) {}

  uint64_t packEggs(std::vector<Egg> eggs, BottomlessBag &bag) { return 0; }

  static void quickSortPar(std::vector<GrainOfSand>::iterator start,
                           std::vector<GrainOfSand>::iterator end,
                           TeamAdventure *team, int shamans) {
    int len = std::distance(start, end);

    if (len >= 2) {
      auto pivot = part(start, end);
      if (shamans > 0 && len > THRESHOLD) {
        auto parallel = (*team).councilOfShamans.enqueue(
            quickSortPar, start, pivot, team, shamans - 2);
        quickSortPar(pivot + 1, end, team, shamans - 2);

        parallel.wait();
      } else {
        quickSortPar(start, pivot, team, 0);
        quickSortPar(pivot + 1, end, team, 0);
      }
    }
  }

  virtual void arrangeSand(std::vector<GrainOfSand> &grains) {
    quickSortPar(grains.begin(), grains.end(), this,
                 static_cast<int>(numberOfShamans));
  }

  virtual Crystal selectBestCrystal(std::vector<Crystal> &crystals) {
    std::vector<std::future<Crystal>> results(numberOfShamans);

    auto begin = crystals.begin();
    uint64_t last = crystals.size() - 1;

    /* Finding Max giving a range for each thread. Corner case for
     * i==numberOfShamans-1 */
    for (uint64_t i = 0; i < numberOfShamans; i++) {
      results[i] =
          councilOfShamans.enqueue(rangeMax, begin + last * i / numberOfShamans,
                                   begin + last * (i + 1) / numberOfShamans -
                                       (i != (numberOfShamans - 1)));
    }

    Crystal max = crystals.front();

    for (uint64_t i = 0; i < numberOfShamans; i++) {
      Crystal res = results[i].get();
      max = (max < res) ? res : max;
    }

    return max;
  }

 private:
  uint64_t numberOfShamans;
  ThreadPool councilOfShamans;

  typedef std::vector<Crystal>::iterator vec_it;

  static Crystal rangeMax(vec_it begin, vec_it end) {
    Crystal max = *begin;

    for (auto iter = begin; iter <= end; ++iter)
      if (max < *iter) max = *iter;

    return max;
  }
};

#endif  // SRC_ADVENTURE_H_
