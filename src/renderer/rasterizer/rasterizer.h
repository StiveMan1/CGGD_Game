#pragma once

#include "resource.h"

#include <functional>
#include <iostream>
#include <linalg.h>
#include <memory>


using namespace linalg::aliases;

namespace cg::renderer
{
	template<typename VB, typename RT>
	class rasterizer
	{
	public:
		rasterizer(){};
		~rasterizer(){};
		void set_render_target(
				std::shared_ptr<resource<RT>> in_render_target,
				std::shared_ptr<resource<float>> in_depth_buffer = nullptr);
		void clear_render_target(
				const RT& in_clear_value, const float in_depth = FLT_MAX);

		void set_vertex_buffer(std::shared_ptr<resource<VB>> in_vertex_buffer);
		void set_index_buffer(std::shared_ptr<resource<unsigned int>> in_index_buffer);

		void set_viewport(size_t in_width, size_t in_height);

		void draw(size_t num_vertexes, size_t vertex_offest);

		std::function<std::pair<float4, VB>(float4 vertex, VB vertex_data)> vertex_shader;
		std::function<cg::color(const VB& vertex_data, const float z)> pixel_shader;

	protected:
		std::shared_ptr<cg::resource<VB>> vertex_buffer;
		std::shared_ptr<cg::resource<unsigned int>> index_buffer;
		std::shared_ptr<cg::resource<RT>> render_target;
		std::shared_ptr<cg::resource<float>> depth_buffer;

		size_t width = 1920;
		size_t height = 1080;

		float edge_function(float2 a, float2 b, float2 c);
		bool depth_test(float z, size_t x, size_t y);
	};

	template<typename VB, typename RT>
	inline void rasterizer<VB, RT>::set_render_target(
			std::shared_ptr<resource<RT>> in_render_target,
			std::shared_ptr<resource<float>> in_depth_buffer)
	{
		if (in_render_target)
			render_target = in_render_target;
		if (in_depth_buffer)
			depth_buffer = in_depth_buffer;
	}

	template<typename VB, typename RT>
	inline void rasterizer<VB, RT>::set_viewport(size_t in_width, size_t in_height)
	{
		width = in_width;
		height = in_height;
	}

	template<typename VB, typename RT>
	inline void rasterizer<VB, RT>::clear_render_target(
			const RT& in_clear_value, const float in_depth)
	{
		if (render_target)
			for (size_t i = 0; i < render_target->get_number_of_elements(); i++)
				render_target->item(i) = in_clear_value;

		if (depth_buffer)
			for (size_t i = 0; i < depth_buffer->get_number_of_elements(); i++)
				depth_buffer->item(i) = in_depth;
	}

	template<typename VB, typename RT>
	inline void rasterizer<VB, RT>::set_vertex_buffer(
			std::shared_ptr<resource<VB>> in_vertex_buffer)
	{
		vertex_buffer = in_vertex_buffer;
	}

	template<typename VB, typename RT>
	inline void rasterizer<VB, RT>::set_index_buffer(
			std::shared_ptr<resource<unsigned int>> in_index_buffer)
	{
		index_buffer = in_index_buffer;
	}

	template<typename VB, typename RT>
	inline void rasterizer<VB, RT>::draw(size_t num_vertexes, size_t vertex_offset)
	{
		std::vector<VB> vertices(3);
		size_t vertex_id = vertex_offset;
		float4 coords;
		float2 vertex_a, vertex_b, vertex_c, min_vertex, bounding_box_begin, max_vertex, bounding_box_end, point;
		float edge, edge0, edge1, edge2, v, u, w, depth;
		size_t u_x, u_y;

		while (vertex_id < vertex_offset + num_vertexes) {
			vertices[0] = vertex_buffer->item(index_buffer->item(vertex_id++));
			vertices[1] = vertex_buffer->item(index_buffer->item(vertex_id++));
			vertices[2] = vertex_buffer->item(index_buffer->item(vertex_id++));

			for (auto& vertex: vertices) {
				coords = {vertex.x, vertex.y, vertex.z, 1.f};
				auto processing_vertex = vertex_shader(coords, vertex);

				vertex.x = processing_vertex.first.x / processing_vertex.first.w;
				vertex.y = processing_vertex.first.y / processing_vertex.first.w;
				vertex.z = processing_vertex.first.z / processing_vertex.first.w;

				vertex.x = (vertex.x + 1.f) * width / 2.f;
				vertex.y = (-vertex.y + 1.f) * height / 2.f;
			}

			{
				vertex_a = float2{vertices[0].x, vertices[0].y};
				vertex_b = float2{vertices[1].x, vertices[1].y};
				vertex_c = float2{vertices[2].x, vertices[2].y};

				min_vertex = min(vertex_a, min(vertex_b, vertex_c));
				bounding_box_begin = round(clamp(min_vertex, float2{0, 0}, float2{(float) (width - 1), (float) (height - 1)}));

				max_vertex = max(vertex_a, max(vertex_b, vertex_c));
				bounding_box_end = round(clamp(max_vertex, float2{0, 0}, float2{(float) (width - 1), (float) (height - 1)}));

				edge = edge_function(vertex_a, vertex_b, vertex_c);
				for (float x = bounding_box_begin.x; x <= bounding_box_end.x; x += 1.f) {
					for (float y = bounding_box_begin.y; y <= bounding_box_end.y; y += 1.f) {
						point = {x, y};
						edge0 = edge_function(vertex_a, vertex_b, point);
						edge1 = edge_function(vertex_b, vertex_c, point);
						edge2 = edge_function(vertex_c, vertex_a, point);

						if (edge0 >= 0.f && edge1 >= 0.f && edge2 >= 0.f) {
							u = edge1 / edge;
							v = edge2 / edge;
							w = edge0 / edge;

							depth = u * vertices[0].z + v * vertices[1].z + w * vertices[2].z;
							u_x = static_cast<size_t>(x);
							u_y = static_cast<size_t>(y);
							if (depth_test(depth, u_x, u_y)) {
								auto pixel_result = pixel_shader(vertices[0], depth);
								render_target->item(u_x, u_y) = RT::from_color(pixel_result);
								if (depth_buffer)
									depth_buffer->item(u_x, u_y) = depth;
							}
						}
					}
				}
			}
		}
	}

	template<typename VB, typename RT>
	inline float
	rasterizer<VB, RT>::edge_function(float2 a, float2 b, float2 c)
	{
		return (c.x - a.x) * (b.y - a.y) - (c.y - a.y) * (b.x - a.x);
	}

	template<typename VB, typename RT>
	inline bool rasterizer<VB, RT>::depth_test(float z, size_t x, size_t y)
	{
		if (!depth_buffer) return true;
		return depth_buffer->item(x, y) > z;
	}

}// namespace cg::renderer