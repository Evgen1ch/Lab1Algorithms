#pragma once

#include <chrono>

class Timer
{
public:
	Timer() : m_duration(0){}
	void start() const{
		m_start = std::chrono::high_resolution_clock::now();
	}
	void end() const{
		m_end = std::chrono::high_resolution_clock::now();
		m_duration = m_end - m_start;
	}
	double getDuration() const
	{
		return m_duration.count();
	}

private:
	mutable std::chrono::time_point<std::chrono::steady_clock> m_start;
	mutable std::chrono::time_point<std::chrono::steady_clock> m_end;
	mutable std::chrono::duration<double> m_duration;
};

