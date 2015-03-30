use inc::Module::Install;

# Define metadata
name            'pAss';
all_from        'lib/PASS/Alignment.pm';
author          'Wesley GOI <wesley@bic.nus.edu.sg>';
perl_version    '5.20';
# Specific dependencies
requires       'Moo'  => '1.007000';
requires       'namespace::clean'  => '0.25';
requires       'autodie' => '2.23';

requires       'Bio::SeqIO' => '1.006922';
requires       'Statistics::Basic' => '1.6611';
requires       'List::MoreUtils' => '0.408';
requires        'experimental' => '0.007';
requires    'Bio::Seq' =>'1.006922';
  #test_requires  'Test::More'  => '0.42';
  #recommends     'Text::CSV_XS'=> '0.50';
  #no_index       'directory'   => 'demos';
  #install_script 'myscript';

WriteAll;
