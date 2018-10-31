=head1 LICENSE
Copyright [2018] EMBL-European Bioinformatics Institute

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

     http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

=cut

package EmailReport;

use strict;
use warnings;
use feature 'say';

use base ('Bio::EnsEMBL::Hive::RunnableDB::NotifyByEmail');
use Data::Dumper;

sub fetch_input {
  my $self = shift;

  my $reports  = $self->param('reports');

  my $subject = "Splice extraction done";
  $self->param('subject', $subject);

  my @text_lines;
  for my $species (sort keys %$reports) {
    my $report = $reports->{$species};
    push @text_lines, $report;
  }
  my $text = join("\n", @text_lines);
  
  $self->param('text', $text);
}

1;
