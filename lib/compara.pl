#!/usr/env perl
use v5.14.00;
use strict;
use warnings;
use Carp;
use autodie qw(:all);
use Readonly;
use Getopt::Long qw(:config no_ignore_case);
use Log::Log4perl qw( :easy ); 
Log::Log4perl->easy_init($WARN); 
my $logger = get_logger(); 

use Bio::EnsEMBL::Registry;
use DBI;

###############################################################################
# MAIN
# Get command line args
my %opt = %{ opt_check() };

Bio::EnsEMBL::Registry->load_all($opt{registry});

# Get the list of genes covered in the db of interest
my $genes1 = get_genes_from_db($opt{db1});
my $genes2 = get_genes_from_db($opt{db2});
$logger->info(sprintf("%s: %d genes", $opt{species1}, scalar @$genes1));
$logger->info(sprintf("%s: %d genes", $opt{species2}, scalar @$genes2));

$genes1 = add_lengths($genes1, $opt{species1});
$genes2 = add_lengths($genes2, $opt{species2});

# Get the correspondence via compara
my $match = get_compara($genes1, $opt{species2});

# Merge the genes and their orthologs
my $merged = merge_matches($genes1, $genes2, $match);
$logger->info("Matched genes: " . scalar(@$merged));

# Print out the result to stdout
print_match($merged);

sub get_genes_from_db {
  my ($inputdb) = @_;

  my @genes;

  # Open SQLite db
  my $dbh = DBI->connect("dbi:SQLite:dbname=$inputdb");

  # Get the list of genes
  my $sth = $dbh->prepare("SELECT * FROM genes");
  $sth->execute();
  for my $gene (values %{ $sth->fetchall_hashref("gene") }) {
    push @genes, $gene;
  }

  return \@genes;
}

sub add_lengths {
  my ($genes, $species) = @_;
  my $tla = Bio::EnsEMBL::Registry->get_adaptor($species, 'core', 'translation');

  $logger->info("Add lengths for $species");
  
  my %lengths;
  my $n = 0;
  for my $tl (@{ $tla->fetch_all() }) {
    $n++; $logger->debug("$n translations...") if $n % 1000 == 0;
    my $length = $tl->length();
    my $gene_id = $tl->stable_id;
    $gene_id =~ s/-P[A-Z]+$//;
    if (not exists $lengths{$gene_id} or $lengths{$gene_id} < $length) {
      $lengths{$gene_id} = $length;
    }
  }

  $logger->debug("Transfer lengths");
  for my $gene (@$genes) {
    $gene->{length} = $lengths{ $gene->{gene} };
  }
  return $genes;
}

sub get_compara {
  my ($genes, $species2) = @_;


  # Get the genomedb for species2
  my $genomedba = Bio::EnsEMBL::Registry->get_adaptor('vb', 'compara', 'GenomeDB');
  my ($genomedb2) = @{ $genomedba->fetch_all_by_name($species2) };
  $logger->info("Recognized $species2 in genomedb table") if $genomedb2;

  my $genea = Bio::EnsEMBL::Registry->get_adaptor('vb', 'compara', 'GeneMember');
  my $homoa = Bio::EnsEMBL::Registry->get_adaptor('vb', 'compara', 'Homology');

  my %matches;
  my $n = 0;
  for my $gene (@$genes) {
    $n++;
    $logger->debug(sprintf("%d/%d", $n, scalar @$genes)) if $n % 1000 == 0;
    #last if $n >= 100;
    my $member = $genea->fetch_by_stable_id($gene->{gene});
    next if not defined $member;
    my $homologies = $homoa->fetch_all_by_Member($member);

    $matches{ $gene->{gene} } = [];
    for my $hom (@$homologies) {
      my $matched = $hom->get_Member_by_GenomeDB($genomedb2);
      if ($matched) {
        for my $match (@$matched) {
          push @{ $matches{ $gene->{gene}} }, $match->stable_id;
        }
      }
    }
  }

  return \%matches;
}

sub merge_matches {
  my ($genes1, $genes2, $matches) = @_;

  my %genes2_hash = map { $_->{gene} => $_ } @$genes2;

  my @merged;
  for my $gene1 (@$genes1) {
    my $matched = $matches->{$gene1->{gene}};
    for my $match (@$matched) {
      $match =~ s/-P[A-Z]+$//;
      my $gene2 = $genes2_hash{$match};
      if (defined $gene2) {
        push @merged, merged($gene1, $gene2);
      } else {
        $logger->warn("No match: $gene1->{gene} vs $match");
      }
    }
  }
  
  return \@merged;
}

sub merged {
  my ($g1, $g2) = @_;
  
  return $g1->{gene} if not $g2;
  my $complete1 = $g1->{introns} > 0 ? $g1->{covered_introns}/$g1->{introns} : 1;
  my $complete2 = $g2->{introns} > 0 ? $g2->{covered_introns}/$g2->{introns} : 2;
  
  my $support = compared($g1, $g2);
  my $introns_diff = compared_introns($g1, $g2);
  my $total_diff = total_diff($support, $introns_diff);
  my @elts = (
    $g1->{gene},
    $g1->{length},
    $g1->{introns},
    $g1->{covered_introns},
    $complete1,
    $g2->{gene},
    $g2->{length},
    $g2->{introns},
    $g2->{covered_introns},
    $complete2,
    $support,
    $introns_diff,
    $total_diff,
  );

  return join("\t", @elts);
}

sub compared {
  my ($g1, $g2) = @_;
  
  my $support1 = support($g1);
  my $support2 = support($g2);

  if ($support1 eq $support2) {
    return "both_" . $support1;
  } else {
    return $support1 . "_vs_" . $support2;
  }
}

sub support {
  my ($g) = @_;

  if ($g->{introns} == 0) {
    return "nointron";
  }
  my $completeness = $g->{covered_introns}/$g->{introns};

  if ($completeness == 0) {
    return "nosupport";
  }
  elsif ($completeness == 1) {
    return "full";
  }
  else {
    return "partial";
  }
}

sub compared_introns {
  my ($g1, $g2) = @_;
  
  my $introns1 = $g1->{introns};
  my $introns2 = $g2->{introns};

  if ($introns1 > $introns2) {
    return "more_introns";
  } elsif ($introns1 < $introns2) {
    return "less_introns";
  } else {
    return "same_introns";
  }
}

sub total_diff {
  my ($support, $introns_diff) = @_;

  if (($support eq 'both_full' or $support eq 'both_nointron')
      and $introns_diff eq 'same_introns') {
    return "perfect_match";
  }
  if ($support eq 'both_full' or $support eq 'both_nointron') {
    return "okay";
  }
  elsif ($support eq 'both_nosupport'
      and $introns_diff eq 'same_introns') {
    return "same_nosupport";
  }
  elsif ($support eq 'both_nosupport') {
    return "bad_support";
  }
  elsif ($support =~ /^full_/) {
    return "better";
  }
  elsif ($support =~ /_full$/) {
    return "worse";
  } else {
    return "other";
  }
}

sub print_match {
  my ($merged) = @_;

  my @header = qw(
    #sp1 len1 intr1 cov1 cov1/intr1
    #sp2 len2 intr2 cov2 cov2/intr2
    support
    introns_diff
    total_diff
  );
  say join("\t", @header);
  for my $line (@$merged) {
    say $line;
  }
}

###############################################################################
# Parameters and usage
sub usage {
  my $error = shift;
  my $help = '';
  if ($error) {
    $help = "[ $error ]\n";
  }
  $help .= <<'EOF';
    Return a list of genes with their match in another species and their splice junction coverage.

    --db1 <path>      : db from which the genes of interest are used
    --db2 <path>      : db compared

    --species1 <str>  : species for db1
    --species2 <str>  : species for db2

    --registry <path> : path to the registry including compara and the species core dbs

    
    --help            : show this help message
    --verbose         : show detailed progress
    --debug           : show even more information (for debugging purposes)
EOF
  print STDERR "$help\n";
  exit(1);
}

sub opt_check {
  my %opt = ();
  GetOptions(\%opt,
    "db1=s",
    "db2=s",
    "species1=s",
    "species2=s",
    "registry=s",
    "help",
    "verbose",
    "debug",
  );

  usage()                if $opt{help};

  usage("db1 and db2 needed") if not $opt{db1} or not $opt{db2};
  usage("species1 and species2 needed") if not $opt{species1} or not $opt{species2};
  usage("Registry file needed") if not $opt{registry};

  Log::Log4perl->easy_init($INFO) if $opt{verbose};
  Log::Log4perl->easy_init($DEBUG) if $opt{debug};
  return \%opt;
}

__END__

