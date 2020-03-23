function random_hairpins = load_random_hairpin_data_20191019(file_random_hairpin)

load(file_random_hairpin);
random_hairpins = [];
random_hairpins.n_bp = n_bp';
random_hairpins.MFE = MFE';
random_hairpins.n_hairpins = n_hairpin';
random_hairpins.loop_size = loop_size';
random_hairpins.fraction_in_stem = fraction_in_stem';
