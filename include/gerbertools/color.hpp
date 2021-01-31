/**
 * MIT License
 *
 * Copyright (c) 2021 Jeroen van Straten
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */

/** \file
 * Stuff related to colors.
 */

#pragma once

namespace gerbertools {

/**
 * Contains stuff related to colors.
 */
namespace color {

/**
 * A color with transparency channel.
 */
struct Color {

    /**
     * Red channel, between 0 and 1.
     */
    float r;

    /**
     * Green channel, between 0 and 1.
     */
    float g;

    /**
     * Blue channel, between 0 and 1.
     */
    float b;

    /**
     * Opacity, between 0 (clear) and 1 (solid).
     */
    float a;

};

/**
 * No color; layer is omitted.
 */
static const Color NONE = {0.0f, 0.0f, 0.0f, 0.0f};

/**
 * Default black color.
 */
static const Color BLACK = {0.0f, 0.0f, 0.0f, 1.0f};

/**
 * Default color for copper.
 */
static const Color COPPER = {0.8f, 0.7f, 0.3f, 1.0f};

/**
 * Default color for HASL/tin surface finish.
 */
static const Color FINISH_TIN = {0.7f, 0.7f, 0.7f, 1.0f};

/**
 * Default color for FR-4 substrate.
 */
static const Color SUBSTRATE = {0.6f, 0.5f, 0.3f, 0.95f};

/**
 * Default color for green soldermask.
 */
static const Color MASK_GREEN = {0.1f, 0.6f, 0.3f, 0.6f};

/**
 * Default color for white soldermask.
 */
static const Color MASK_WHITE = {0.9f, 0.93f, 1.00f, 0.8f};

/**
 * Default color for white silkscreen.
 */
static const Color SILK_WHITE = {0.9f, 0.9f, 0.9f, 0.9f};

/**
 * Default color for black silkscreen.
 */
static const Color SILK_BLACK = {0.1f, 0.1f, 0.1f, 0.9f};

} // namespace color
} // namespace gerbertools
