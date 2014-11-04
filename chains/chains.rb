#! /usr/bin/ruby

require "matrix"

class Chain

	attr_reader :lattice_sites
	attr_reader :ions
	attr_reader :all_sites
	attr_reader :displacement
	attr_reader :timesteps

	def initialize( timestep, ion, sites, displacement )
		@timesteps = [ timestep ]
		@ions = [ ion ]
		@all_sites = sites
		@lattice_sites = [ sites[0], sites[-1] ]
		@displacement = Vector[*displacement]
	end

	def extended_by?( new_chain )
		return @lattice_sites[-1] == new_chain.lattice_sites[0]
	end

	def append( new_chain )
		@ions = @ions | new_chain.ions
		@timesteps += new_chain.timesteps
		@lattice_sites += new_chain.lattice_sites[1..-1]
		@all_sites += new_chain.all_sites[1..-1]
		@displacement += new_chain.displacement
	end

	def closed_loop?
		@lattice_sites[0] == @lattice_sites[-1]
	end

	def size
		closed_loop? ? @lattice_sites.size - 1 : @lattice_sites.size
	end

	def output # e.g. => "3: [c] [ 0.0] -265 -484 -260 -265: #109 #204 #459"
		"#{"%3d" % size}: [#{closed_loop? ? 'c' : 'o'}] [#{"%4.1f" % dr}] #{@lattice_sites.join(' ')}: ##{@ions.sort.join(' #')} "
	end

	def <=>( other )
		size <=> other.size
	end

	def dr
		Math.sqrt( @displacement.to_a.inject(0) { |result, element| result += element * element } )
	end

	def duration
		@timesteps.max - @timesteps.min
	end

end

diff_events = File.open( "diff_events.list", 'r' )

def extract_data_from( line ) # line has format of e.g. => "ion  440:      290 =>   427 ( 738 -399 -439 ) [ -0.00  5.12 -7.09 ]"
	data = line.split(/[:\s]+/)
	timestep = data[0].to_i
	ion = data[2].to_i
	sites = data[7..-7].unshift( data[3] ).push( data[5] ).map{ |s| s.to_i }
	displacement = data[-4..-2].map{ |s| s.to_f }
	return timestep, ion, sites, displacement
end

def find_overlap_between( chains )
	finished = false
	until finished
		finished = true
		chains.each do |chain_i|
			chains.each do |chain_j|
				next if chain_i == chain_j
				if chain_i.extended_by?( chain_j )
					chain_i.append( chain_j )
					chains.delete( chain_j )
					finished = false
				elsif chain_j.extended_by?( chain_i )
					chain_j.append( chain_i )
					chains.delete( chain_i )
					finished = false
					break
				end
			end
		end
	end
	return chains
end

def distribution_of( chains )
  distribution = {}
  chains.each do |chain|
    if distribution[chain.size].nil?
      distribution[chain.size] = 1
    else
      distribution[chain.size] += 1
    end
  end    
  return distribution
end

chains = diff_events.readlines.collect{ |line| Chain.new( *extract_data_from( line ) ) }
find_overlap_between( chains )
chains.sort.each{ |chain| puts "(#{chain.duration}) #{chain.output}" }

chains_out = File.new( "chains.out", 'w' )
distribution_of( chains ).sort.each{ |size,number| chains_out.puts "#{size} #{number}"}