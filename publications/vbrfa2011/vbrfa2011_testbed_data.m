function vbrfa2011_testbed_data()

rand('state', 10);
randn('state', 10);

data = testbed_preprocess(testbed_loaddata());

[M,N] = size(data.observations);
Itrain = rand(M,N) > 0.2;
Itest = ~Itrain;

data_train = data;
data_train.observations(~Itrain) = NaN;

data_test = data;
data_test.observations(~Itest) = NaN;

folder = '/share/climate/jluttine/testbed/';
file_train = [folder, 'testbed_vbrfa2011_traindata'];
file_test = [folder, 'testbed_vbrfa2011_testdata'];
save(file_train, '-struct', 'data_train');
save(file_test, '-struct', 'data_test');

fprintf('Saved train and test data into %s and %s\n', file_train, file_test);
