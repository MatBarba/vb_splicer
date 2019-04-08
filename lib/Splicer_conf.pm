package Splicer_conf;

use strict;
use warnings;

use base ('Bio::EnsEMBL::Hive::PipeConfig::HiveGeneric_conf');

use Bio::EnsEMBL::Hive::PipeConfig::HiveGeneric_conf;
use Bio::EnsEMBL::Hive::Version 2.4;

use File::Spec::Functions qw(catdir);
use FindBin;

sub pipeline_wide_parameters {
    my ($self) = @_;

    return {
        %{$self->SUPER::pipeline_wide_parameters},
        'debug' => $self->o('debug'),
        'tmp_dir' => $self->o('tmp_dir'),
    };
}

sub default_options {
  my ($self) = @_;
  return {
    %{ $self->SUPER::default_options() },

    pipeline_name => 'vb_splicer',
    email => $ENV{USER} . '@ebi.ac.uk',

    do_donor_acceptor => 0,
    rest_server => undef,
    do_paralogs => 0,
    homologs_dir => undef,
    paralogs_dir => undef,
    
    # Species factory
    species => [],
    antispecies => [],
    division => [],
    run_all => 0,
    meta_filters => {},
    
    force_gtf => 0,
    force_splice_db => 0,
    force_merge => 0,
    do_not_retag => 0,
    coverage => 1,

    debug => 0,
  };
}

sub hive_meta_table {
  my ($self) = @_;
  return {
    %{$self->SUPER::hive_meta_table},       # here we inherit anything from the base class
    'hive_use_param_stack'  => 1,           # switch on the new param_stack mechanism
  };
}

sub pipeline_analyses {
  my ($self) = @_;

  return
  [
    {
      -logic_name => 'Start',
      -module     => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
      -input_ids  => [{}],
      -flow_into  => {
        '1->A' => 'Para_parser',
        'A->1' => 'Email_report',
      },
      -rc_name    => 'normal',
      -meadow_type       => 'LSF',
    },

    {
      -logic_name => 'Email_report',
      -module            => 'EmailReport',
      -parameters => {
        email => $self->o('email'),
      },
      -rc_name    => 'normal',
      -meadow_type       => 'LSF',
    },

    {
      -logic_name => 'Para_parser',
      -module            => 'ParaParser',
      -language   => 'python3',
      -parameters => {
        do_paralogs      => $self->o('do_paralogs'),
        homologs_dir      => $self->o('homologs_dir'),
        paralogs_dir      => $self->o('paralogs_dir'),
      },
      -flow_into  => 'Species_factory',
      -max_retry_count => 0,
      -rc_name    => 'normal',
      -meadow_type       => 'LSF',
    },

    {
      -logic_name => 'Species_factory',
      -module     => 'Bio::EnsEMBL::Production::Pipeline::Common::SpeciesFactory',
      -parameters => {
        species      => $self->o('species'),
        antispecies  => $self->o('antispecies'),
        run_all      => $self->o('run_all'),
        meta_filters => {},
      },
      -max_retry_count => 0,
      -flow_into  => {
        '2' => 'SpeciesFilter',
      },
      -rc_name    => 'normal',
      -meadow_type       => 'LSF',
    },

    {
      -logic_name => 'SpeciesFilter',
      -module            => 'SpeciesFilter',
      -language   => 'python3',
      -parameters => {
        bam_dir => $self->o('bam_dir'),
      },
      -flow_into  => {
        '2->A' => 'Files',
        'A->2' => 'Report',
      },
      -rc_name    => 'normal',
      -meadow_type       => 'LSF',
    },

    {
      -logic_name => 'Report',
      -module            => 'Report',
      -language   => 'python3',
      -flow_into  => {
        '2' => '?accu_name=reports&accu_address={species}&accu_input_variable=report',
      },
      -analysis_capacity => 1,
      -rc_name    => 'normal',
      -meadow_type       => 'LSF',
    },

    {
      -logic_name => 'Files',
      -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
      -parameters        => {
        tmp_dir => $self->o('tmp_dir'),
        cmd => "mkdir -p #tmp_dir#",
      },
      -flow_into  => {
        '1' => ['GTF_factory', 'Chrom_sizes'],
      },
      -rc_name    => 'normal',
      -meadow_type       => 'LSF',
    },

    {
      -logic_name        => 'Chrom_sizes',
      #-module            => 'ChromSizes',
      -module     => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
      #-language   => 'python3',
      -parameters        => {
        tmp_dir => $self->o('tmp_dir'),
        rest_server => $self->o('rest_server'),
        sizes => "#tmp_dir#/#species#.sizes",
      },
      -flow_into  => {
        '1' => {'BigWig_merge' => { sizes => "#sizes#" } },
      },
      -analysis_capacity => 1,
      -max_retry_count   => 0,
      -failed_job_tolerance => 0,
      -rc_name    => 'normal',
      -meadow_type       => 'LSF',
    },


    {
      -logic_name        => 'BigWig_merge',
      -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
      -parameters        => {
        bigwig_dir => $self->o('bigwig_dir'),
        bigwig_list => "#bigwig_dir#/#species#/*",
        output => "#tmp_dir#/#species#.bed",
        cmd => "bigWigMerge #bigwig_list# #output#",
      },
      -flow_into  => {
        '1' => {'Bed_to_BigWig' => { bed => "#output#", sizes => "#sizes#" } },
      },
      -analysis_capacity => 0,
      -max_retry_count   => 0,
      -failed_job_tolerance => 0,
      -rc_name    => 'normal',
      -meadow_type       => 'LSF',
    },

    {
      -logic_name        => 'Bed_to_BigWig',
      -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
      -parameters        => {
        file => "#tmp_dir#/#species#.bed",
        cmd => "bedGraphToBigWig #bed# #sizes# #file#",
      },
      -flow_into  => {
        '1' => '?accu_name=files&accu_address={species}[]&accu_input_variable=file',
      },
      -analysis_capacity => 0,
      -max_retry_count   => 0,
      -failed_job_tolerance => 0,
      -rc_name    => 'normal',
      -meadow_type       => 'LSF',
    },

    {
      -logic_name        => 'GTF_factory',
      -module            => 'GtfFactory',
      -language   => 'python3',
      -parameters        => {
        gtf_dir       => $self->o('gtf_dir'),
        out_file_stem => 'gtf',
        force_gtf     => $self->o('force_gtf'),
      },
      -analysis_capacity => 0,
      -max_retry_count   => 0,
      -failed_job_tolerance => 0,
      -flow_into  => {
        '2->A' => 'GTF_dumper',
        'A->3' => 'SpliceDB_factory',
      },
      -rc_name    => 'normal',
      -meadow_type       => 'LSF',
    },

    {
      -logic_name        => 'GTF_dumper',
      -module            => 'Bio::EnsEMBL::EGPipeline::FileDump::GTFDumper',
      -max_retry_count   => 0,
      -parameters        => {
        data_type   => 'basefeatures',
        db_type     => 'core',
        results_dir => $self->o('gtf_dir'),
        out_file_stem => 'gtf',
      },
      -analysis_capacity => 20,
      -rc_name           => 'bigmem',
      -meadow_type       => 'LSF',
    },
    
    {
      -logic_name        => 'SpliceDB_factory',
      -module            => 'SpliceDBFactory',
      -language   => 'python3',
      -parameters        => {
        bam_dir     => $self->o('bam_dir'),
        splice_dir     => $self->o('splice_dir'),
        force_splice_db => $self->o('force_splice_db'),
      },
      -max_retry_count   => 0,
      -flow_into  => {
        '2->A' => 'Extract_splices',
        'A->3' => 'Merge_splices',
        '3' => 'Merge_splices',
      },
      -rc_name    => 'normal',
      -meadow_type       => 'LSF',
    },

    {
      -logic_name => 'Extract_splices',
      -module     => 'ExtractSplices',
      -language   => 'python3',
      
      -analysis_capacity => 30,
      -max_retry_count => 0,
      -meadow_type       => 'LSF',
      -rc_name    => 'bigmem',
      #-flow_into  => {
      #  '1' => '?accu_name=splice_dbs&accu_address={species}[]&accu_input_variable=splice_db',
      #}   
    },
 
    {
      -logic_name => 'Merge_splices',
      -module     => 'MergeSplices',
      -language   => 'python3',
      -parameters        => {
        splice_dir     => $self->o('splice_dir'),
        force_merge => $self->o('force_merge'),
      },
      -analysis_capacity => 10,
      -max_retry_count => 0,
      -meadow_type       => 'LSF',
      -rc_name    => 'bigmem',
      -flow_into  => {
        '-1' => 'Merge_splices_highmem',
        '2' => 'Tagger',
      }
    },
 
    {
      -logic_name => 'Merge_splices_highmem',
      -module     => 'MergeSplices',
      -language   => 'python3',
      -parameters        => {
        splice_dir     => $self->o('splice_dir'),
        force_merge => $self->o('force_merge'),
      },
      -analysis_capacity => 10,
      -max_retry_count => 0,
      -meadow_type       => 'LSF',
      -rc_name    => 'biggermem',
      -flow_into  => {
        '2' => 'Tagger_highmem',
      }
    },

    {
      -logic_name => 'Tagger',
      -module     => 'Tagger',
      -language   => 'python3',
      -parameters        => {
        paralogs_dir =>  $self->o('paralogs_dir'),
        do_not_retag => $self->o('do_not_retag'),
        rest_server => $self->o('rest_server'),
        do_donor_acceptor => $self->o('do_donor_acceptor'),
      },
      -analysis_capacity => 20,
      -max_retry_count => 0,
      -rc_name    => 'bigmem',
      -meadow_type       => 'LSF',
      -flow_into  => {
        '1' => 'Output_factory',
        '-1' => 'Tagger_highmem',
      }
    },

    {
      -logic_name => 'Tagger_highmem',
      -module     => 'Tagger',
      -language   => 'python3',
      -parameters        => {
        paralogs_dir =>  $self->o('paralogs_dir'),
        do_not_retag => $self->o('do_not_retag'),
        rest_server => $self->o('rest_server'),
        do_donor_acceptor => $self->o('do_donor_acceptor'),
      },
      -analysis_capacity => 20,
      -max_retry_count => 0,
      -rc_name    => 'biggermem',
      -meadow_type       => 'LSF',
      -flow_into  => {
        '1' => 'Output_factory',
      }
    },

    {
      -logic_name => 'Output_factory',
      -module     => 'OutputFactory',
      -language   => 'python3',
      -parameters        => {
        summary_dir     => $self->o('summary_dir'),
        ftp_summary_dir     => $self->o('ftp_summary_dir'),
        json_dir    => $self->o('json_dir'),
        rest_server => $self->o('rest_server'),
      },
      -analysis_capacity => 5,
      -max_retry_count => 0,
      -flow_into  => {
        '2' => 'Create_GFF',
        '3' => 'Create_Apollo_json',
      },
      -rc_name    => 'normal',
      -meadow_type       => 'LSF',
    },

    {
      -logic_name => 'Create_GFF',
      -module     => 'CreateGFF',
      -language   => 'python3',
      -parameters        => {
        category    => '#category#',
        coverage    => '#coverage#',
        track_file  => '#track_file#',
      },
      -hive_capacity => 20,
      -max_retry_count => 0,
      -rc_name    => 'bigmem',
      -meadow_type       => 'LSF',
      -flow_into  => {
        '2' => '?accu_name=files&accu_address={species}[]&accu_input_variable=file',
        '-1' => 'Create_GFF_highmem',
      }
    },

    {
      -logic_name => 'Create_GFF_highmem',
      -module     => 'CreateGFF',
      -language   => 'python3',
      -parameters        => {
        category    => '#category#',
        coverage    => '#coverage#',
        track_file  => '#track_file#',
      },
      -hive_capacity => 20,
      -max_retry_count => 0,
      -rc_name    => 'biggermem',
      -meadow_type       => 'LSF',
      -flow_into  => {
        '2' => '?accu_name=files&accu_address={species}[]&accu_input_variable=file',
      }
    },

    {
      -logic_name => 'Create_Apollo_json',
      -module     => 'Apollo_track',
      -language   => 'python3',
      -parameters        => {
        label        => '#label#',
        description  => '#description#',
        track_file  => '#track_file#',
        json_file  => '#json_file#',
        version  => '#version#',
      },
      -analysis_capacity => 1,
      -max_retry_count => 0,
      -rc_name    => 'normal',
      -meadow_type       => 'LSF',
      -flow_into  => {
        #  '2' => '?accu_name=gffs&accu_address={species}[]&accu_input_variable=gff',
      }
    },
  ];
}

sub resource_classes {

  my ($self) = @_;
  
  my $reg_requirement = '--reg_conf ' . $self->o('reg_conf');

  return {
    'default'           => {
      'LOCAL' => ['', $reg_requirement],
    },
    'normal'            => {'LSF' => ['-q production-rh7 -M  1000 -R "rusage[mem=1000]"', $reg_requirement]},
    'bigmem'           => {'LSF' => ['-q production-rh7 -M  4000 -R "rusage[mem=4000]"', $reg_requirement]},
    'biggermem'           => {'LSF' => ['-q production-rh7 -M  32000 -R "rusage[mem=32000]"', $reg_requirement]},
  }
}

1;
