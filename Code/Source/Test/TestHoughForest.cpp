#include "../Common.hpp"
#include "../Algorithms/HoughForest.hpp"
#include "../Math.hpp"
#include "../Matrix.hpp"
#include "../UnorderedMap.hpp"
#include <boost/algorithm/string.hpp>
#include <fstream>
#include <sstream>

using namespace std;
using namespace Thea;
using namespace Algorithms;

bool testHoughForestFile(string const & path);

int
main(int argc, char * argv[])
{
  if (argc < 2)
  {
    cout << "Usage: " << argv[0] << " <training-file>" << endl;
    return 0;
  }

  try
  {
    if (!testHoughForestFile(argv[1]))
      return -1;
  }
  catch(...)
  {
    cerr << "An error occurred" << endl;
    return -1;
  }

  return 0;
}

template <typename ArrayT>
string
arrayToString(ArrayT const & arr)
{
  ostringstream oss; oss << '[';
  for (long i = 0; i < (long)arr.size(); ++i)
  {
    if (i > 0) oss << ", ";
    oss << arr[i];
  }
  oss << ']';

  return oss.str();
}

template <typename IteratorT>
string
rangeToString(IteratorT begin, IteratorT end)
{
  ostringstream oss; oss << '[';
  for (IteratorT ii = begin; ii != end; ++ii)
  {
    if (ii != begin) oss << ", ";
    oss << *ii;
  }
  oss << ']';

  return oss.str();
}

/**
 * Training data for Hough forest. Consists of a list of points, each with:
 *   (a) class ID,
 *   (b) Hough vote parameters, and
 *   (c) feature vector.
 *
 * All points are assumed to have the same number of Hough parameters. The input class ID's are mapped to a contiguous set of
 * Hough class ID's in the range [0, num_classes - 1]. The inputClassToHoughClass() function can be used to apply this mapping,
 * and houghClassToInputClass() to invert it. 0 is always the background class (in both the input and Hough cases).
 */
class ExampleSet: public HoughForest::TrainingData
{
  private:
    typedef TheaUnorderedMap<long, long> IDMap;

  public:
    THEA_DEF_POINTER_TYPES(ExampleSet, Thea::shared_ptr, Thea::weak_ptr)

    /** Constructor. */
    ExampleSet(long num_classes_ = 0) : num_classes(num_classes_) {}

    /** Construct with known parameters. */
    ExampleSet(long num_classes_, long num_features_, long num_vote_params_, long num_examples_)
    {
      resize(num_classes_, num_features_, num_vote_params_, num_examples_);
    }

    /** Set the data for a single example point. */
    void setExample(long example_index, double const * example_features, long example_class, double const * example_self_vote)
    {
      alwaysAssertM(example_index >= 0 && example_index < (long)classes.size(), "Example index out of bounds");

      features.setColumn(example_index, example_features);
      classes[(size_t)example_index] = example_class;
      votes.setColumn(example_index, example_self_vote);
    }

    /** Convert an input class ID (read from file) to a Hough class ID in the range [0, num_classes - 1]. */
    long inputClassToHoughClass(long c) const
    {
      IDMap::const_iterator existing = input_class_to_hough_class.find(c);
      if (existing != input_class_to_hough_class.end())
        return existing->second;
      else
        return -1;
    }

    /** Convert a Hough class ID in the range [0, num_classes - 1] to an input class ID (read from file). */
    long houghClassToInputClass(long c) const
    {
      IDMap::const_iterator existing = hough_class_to_input_class.find(c);
      if (existing != hough_class_to_input_class.end())
        return existing->second;
      else
        return -1;
    }

    /**
     * Load the set of example points from a disk file. The file is assumed to have one example point per line, each with the
     * following structure:
     *
     * [model_id] [class_id] [cx cy cz] [features...],
     *
     * where [cx cy cz] are the 3 Hough vote parameters, i.e. the offset to the object center.
     */
    bool load(std::string const & path)
    {
      THEA_CONSOLE << "Reading example set from: " << path;

      long num_vote_params = 3;

      // Read number of examples, features
      long num_examples = 0, num_features = 0;
      hough_class_to_input_class.clear();
      input_class_to_hough_class.clear();
      {
        std::ifstream in(path.c_str());
        if (!in)
        {
          THEA_ERROR << "Could not read file: " << path;
          return false;
        }

        long next_non_background_class_id = 1;

        std::string line;
        while (std::getline(in, line))
        {
          boost::trim(line);
          if (!line.empty())
          {
            num_examples++;

            std::istringstream line_in(line);

            long class_id, model_id;
            line_in >> class_id >> model_id;

            if (inputClassToHoughClass(class_id) < 0)
            {
              long hough_class_id = 0;
              if (class_id != 0)  // foreground class
                hough_class_id = next_non_background_class_id++;

              input_class_to_hough_class[class_id] = hough_class_id;
              hough_class_to_input_class[hough_class_id] = class_id;
            }

            if (num_examples == 1)  // first example
            {
              double x;
              num_features = -num_vote_params;  // ignore fields for vote parameters
              while (line_in >> x)
                num_features++;
            }
          }
        }
      }

      THEA_CONSOLE << "Example set '" << path << '\'';
      THEA_CONSOLE << "  #classes     : " << hough_class_to_input_class.size();
      THEA_CONSOLE << "  #features    : " << num_features;
      THEA_CONSOLE << "  #vote-params : " << num_vote_params;
      THEA_CONSOLE << "  #examples    : " << num_examples;

      // Read actual data
      resize((long)hough_class_to_input_class.size(), num_features, num_vote_params, num_examples);
      {
        std::ifstream in(path.c_str());
        TheaArray<double> example_features((size_t)num_features);
        TheaArray<double> example_self_vote((size_t)num_vote_params);
        long class_id, model_id;
        for (long i = 0; i < num_examples; ++i)
        {
          in >> class_id >> model_id;

          for (long j = 0; j < num_vote_params; ++j)
            in >> example_self_vote[(size_t)j];

          for (long j = 0; j < num_features; ++j)
            in >> example_features[(size_t)j];

          setExample(i, &example_features[0], inputClassToHoughClass(class_id), &example_self_vote[0]);
        }
      }

      THEA_CONSOLE << "Finished reading example set from: " << path;

      return true;
    }

    /** Get the number of example points. */
    long numExamples() const { return (long)classes.size(); }

    /** Get the number of classes (including the background class, always assumed to have ID 0. */
    long numClasses() const { return (long)input_class_to_hough_class.size(); }

    /** Get the number of features. */
    long numFeatures() const { return features.numRows(); }

    /** Get the number of vote parameters for a given class. */
    long numVoteParameters(long class_index) const { return votes.numRows(); }

    /** Get the value of a particular feature for all examples. */
    void getFeatures(long feature_index, double * values) const
    {
      double const * row_start = &features(feature_index, 0);
      std::memcpy(values, row_start, (size_t)features.numColumns() * sizeof(double));
    }

    /** Get the value of a particular feature for a selected set of examples. */
    void getFeatures(long feature_index, long num_selected_examples, long const * selected_examples, double * values) const
    {
      double const * row_start = &features(feature_index, 0);
      for (long i = 0; i < num_selected_examples; ++i)
        values[i] = row_start[selected_examples[i]];
    }

    /** Get the feature vector of a single example. */
    void getExampleFeatures(long example_index, double * values) const
    {
      for (long i = 0; i < features.numRows(); ++i)
        values[i] = features(i, example_index);
    }

    /** Get the Hough class ID's (in the range [0, num_classes - 1]) of all examples. */
    void getClasses(long * classes_) const
    {
      std::memcpy(classes_, &classes[0], classes.size() * sizeof(double));
    }

    /** Get the Hough class ID's (in the range [0, num_classes - 1]) of a selected set of examples. */
    void getClasses(long num_selected_examples, long const * selected_examples, long * classes_) const
    {
      for (long i = 0; i < num_selected_examples; ++i)
        classes_[i] = classes[(size_t)selected_examples[i]];
    }

    /** Get the Hough class ID (in the range [0, num_classes - 1]) of a single example. */
    long getClass(long example_index) const
    {
      return classes[(size_t)example_index];
    }

    /** Get the Hough vote for a single example. */
    void getSelfVote(long example_index, double * params) const
    {
      double const * col_start = &votes(0, example_index);
      std::memcpy(params, col_start, (size_t)votes.numRows() * sizeof(double));
    }

  private:
    /** Resize the set according to a given set of parameters. */
    void resize(long num_classes_, long num_features_, long num_vote_params_, long num_examples_)
    {
      num_classes = num_classes_;
      classes.resize((size_t)num_examples_);
      features.resize(num_features_, num_examples_);
      votes.resize(num_vote_params_, num_examples_);
    }

    long num_classes;
    Matrix<double, MatrixLayout::ROW_MAJOR> features;
    Matrix<double, MatrixLayout::COLUMN_MAJOR> votes;
    TheaArray<long> classes;
    IDMap input_class_to_hough_class;
    IDMap hough_class_to_input_class;

}; // class ExampleSet

/** Called when a vote is cast. */
class VoteCallback : public HoughForest::VoteCallback
{
  public:
    VoteCallback(ExampleSet const * training_data_) : training_data(training_data_) {}

    void operator()(Vote const & vote)
    {
      THEA_CONSOLE << "Hough vote:";
      THEA_CONSOLE << "  - target class : " << training_data->houghClassToInputClass(vote.getTargetClassID());
      THEA_CONSOLE << "  - vote params  : " << rangeToString(vote.getParameters(), vote.getParameters() + vote.numParameters());
      THEA_CONSOLE << "  - weight       : " << vote.getWeight();
    }

  private:
    ExampleSet const * training_data;

}; // class VoteCallback

bool
testHoughForestFile(string const & path)
{
  static int const NUM_TREES = 3;
  static int const NUM_POINTS_PER_CLASS = 3;
  static int const NUM_VOTES_PER_POINT = 3;

  THEA_CONSOLE << "======================================";
  THEA_CONSOLE << "Reading features from file";
  THEA_CONSOLE << "======================================\n";

  ExampleSet training_data;
  if (!training_data.load(path))
    return false;

  THEA_CONSOLE << "\n======================================";
  THEA_CONSOLE << "Training Hough forest";
  THEA_CONSOLE << "======================================\n";

  HoughForest::Options opts = HoughForest::Options::defaults();
  opts.setVerbose(1);  // 2/3 for more debugging output, 0 to disable
  // opts.setMaxDepth(100);
  opts.setMaxLeafElements(3);
  // opts.setMaxDominantFraction(1);
  // opts.setMaxCandidateFeatures(100000000);
  opts.setNumFeatureExpansions(3);
  // opts.setMaxCandidateThresholds(10);
  opts.setProbabilisticSampling(true);

  TheaArray<long> num_vote_params((size_t)training_data.numClasses());
  for (size_t i = 0; i < num_vote_params.size(); ++i)
    num_vote_params[i] = training_data.numVoteParameters((long)i);

  HoughForest hf(training_data.numClasses(), training_data.numFeatures(), &num_vote_params[0], opts);
  hf.train(NUM_TREES, training_data);

  hf.dumpToConsole();

  THEA_CONSOLE << "\n======================================";
  THEA_CONSOLE << "Testing Hough forest";
  THEA_CONSOLE << "======================================";

  // Change the debugging level here if necessary
  hf.setVerbose(1);

  VoteCallback callback(&training_data);
  TheaArray<double> training_vote, test_vote;
  TheaArray<double> features((size_t)training_data.numFeatures());

  for (long c = 1; c < training_data.numClasses(); ++c)
  {
    long inpc = training_data.houghClassToInputClass(c);

    THEA_CONSOLE << "\n======================================";
    THEA_CONSOLE << "Voting for class " << inpc;
    THEA_CONSOLE << "======================================";

    for (long i = 0; i < NUM_POINTS_PER_CLASS; ++i)
    {
      long index = std::rand() % training_data.numExamples();
      long num_vote_params = training_data.numVoteParameters(c);

      training_vote.resize((size_t)num_vote_params);
      test_vote.resize((size_t)num_vote_params);

      training_data.getSelfVote(index, &training_vote[0]);
      training_data.getExampleFeatures(index, &features[0]);

      THEA_CONSOLE << "\n################################################################";
      THEA_CONSOLE << "Testing training example " << index;
      THEA_CONSOLE << "  - features    : " << arrayToString(features);
      THEA_CONSOLE << "  - vote params : " << arrayToString(training_vote);
      THEA_CONSOLE <<   "################################################################\n";

      hf.voteSelf(c, &features[0], NUM_VOTES_PER_POINT, callback);
    }
  }

  return true;
}
