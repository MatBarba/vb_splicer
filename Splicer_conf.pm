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
    };
}

sub default_options {
  my ($self) = @_;
  return {
    %{ $self->SUPER::default_options() },

    pipeline_name => 'vb_splicer',
    email => $ENV{USER} . '@ebi.ac.uk',
    
    # Species factory
    species => [],
    antispecies => [],
    division => [],
    run_all => 0,
    meta_filters => {},
    
    force_gtf => 0,

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
        '1->A' => 'Species_factory',
        'A->1' => 'Email_report',
      },
      -rc_name    => 'default',
      -meadow_type       => 'LOCAL',
    },

    {
      -logic_name => 'Email_report',
      -module            => 'EmailReport',
      -parameters => {
        email => $self->o('email'),
      },
      -rc_name    => 'default',
      -meadow_type       => 'LOCAL',
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
      -meadow_type       => 'LOCAL',
      -max_retry_count => 0,
      -rc_name    => 'default',
      -flow_into  => {
        '2' => 'SpeciesFilter',
      },
    },

    {
      -logic_name => 'SpeciesFilter',
      -module            => 'SpeciesFilter',
      -language   => 'python3',
      -parameters => {
        bam_dir => $self->o('bam_dir'),
      },
      -rc_name    => 'default',
      -meadow_type       => 'LOCAL',
      -flow_into  => {
        '2->A' => 'GTF_factory',
        'A->2' => 'Report',
      },
    },

    {
      -logic_name => 'Report',
      -module            => 'Report',
      -language   => 'python3',
      -rc_name    => 'default',
      -meadow_type       => 'LOCAL',
      -flow_into  => {
        '2' => '?accu_name=reports&accu_address={species}&accu_input_variable=report',
      },
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
      -max_retry_count   => 0,
      -failed_job_tolerance => 0,
      -meadow_type       => 'LOCAL',
      -rc_name           => 'default',
      -flow_into  => {
        '2->A' => 'GTF_dumper',
        'A->3' => 'BAM_factory',
      },
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
      -logic_name        => 'BAM_factory',
      -module            => 'BamFactory',
      -language   => 'python3',
      -parameters        => {
        bam_dir     => $self->o('bam_dir'),
      },
      -max_retry_count   => 0,
      -meadow_type       => 'LOCAL',
      -rc_name           => 'default',
      -flow_into  => {
        '2->A' => 'Extract_splices',
        'A->1' => 'Merge_splices',
      },
    },

    {
      -logic_name => 'Extract_splices',
      -module     => 'ExtractSplices',
      -language   => 'python3',
      -parameters        => {
        splice_dir     => $self->o('splice_dir'),
      },
      
      -analysis_capacity => 1,
      -max_retry_count => 0,
      -meadow_type       => 'LSF',
      -rc_name    => 'bigmem',
      -flow_into  => {
        '2' => '?accu_name=splice_dbs&accu_address={species}[]&accu_input_variable=splice_db',
      }   
    },
 
    {
      -logic_name => 'Merge_splices',
      -module     => 'MergeSplices',
      -language   => 'python3',
      -parameters        => {
        splice_dir     => $self->o('splice_dir'),
      },
      -analysis_capacity => 1,
      -max_retry_count => 0,
      -meadow_type       => 'LSF',
      -rc_name    => 'bigmem',
      -flow_into  => {
        '2' => 'Create_GFF',
      }
    },
    
    {
      -logic_name => 'Create_GFF',
      -module     => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
      -rc_name    => 'default',
      -meadow_type       => 'LOCAL',
      -analysis_capacity => 1,
    },
# 
#     {
#       -logic_name => 'Create_GFF',
#       -module     => 'CreateGFF',
#       -language   => 'python3',
#       -max_retry_count => 0,
#       -rc_name    => 'default',
#       -meadow_type       => 'LOCAL',
#       -flow_into  => {
#         '2' => '?accu_name=gffs&accu_address={species}[]&accu_input_variable=gff',
#       }
#     },
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
  }
}

1;
